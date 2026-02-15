import sys
import time
import numpy as np
import matplotlib.pyplot as plt

# --- COLOR CODES ---
RED = "\033[91m"
RESET = "\033[0m"

# --- SIMULATION PARAMETERS ---
SYSTEM_FREQ = 50.0   # Nominal Frequency (Hz)
H_CONST = 23.64      # System Inertia Constant [cite: 10]
TIME_STEP = 0.5      # Simulation step (seconds)
TRIP_TIME = 5.0      # Time of generator trip
TOTAL_TIME = 20.0    # Total simulation duration

# =================================================================
# 1. SYSTEM DATA & Y-BUS
# =================================================================
def get_ieee_9bus_data():
    # IEEE 9-Bus System Data [cite: 2, 17]
    # format: [id, type(1=Slack, 2=PV, 3=PQ), V, angle_deg, Pg, Qg, Pl, Ql, P_max]
    bus_configs = [
        [1, 1, 1.040, 0.0,   0.764, 0.275, 0.0,   0.0,   1.00], # Slack (P_max limited to 1.0pu)
        [2, 2, 1.025, 8.89,  1.630, 0.073, 0.0,   0.0,   3.00], # Gen 2
        [3, 2, 1.025, 4.24,  0.850, -0.098, 0.0,  0.0,   1.50], # Gen 3
        [4, 3, 1.026, 27.64, 0.0,   0.0,   0.0,   0.0,   0.0],  
        [5, 3, 0.996, 25.81, 0.0,   0.0,   1.25,  0.50,  0.0],  # Load A [cite: 19]
        [6, 3, 1.012, 26.01, 0.0,   0.0,   0.90,  0.30,  0.0],  # Load B
        [7, 3, 1.025, 33.33, 0.0,   0.0,   0.0,   0.0,   0.0],
        [8, 3, 1.015, 30.25, 0.0,   0.0,   1.00,  0.35,  0.0],  # Load C
        [9, 3, 1.032, 31.55, 0.0,   0.0,   0.0,   0.0,   0.0]
    ]
    
    bus_data = []
    for b in bus_configs:
        bus_data.append({
            'id': b[0], 'type': b[1], 'V': b[2], 'theta': np.radians(b[3]),
            'Pg': b[4], 'Pl': b[6], 'Qg': b[5], 'Ql': b[7],
            'P_spec': b[4] - b[6], 'Q_spec': b[5] - b[7], 'P_max': b[8]
        })

    line_configs = [
        [1, 4, 0.0000, 0.0576, 0.0000], [4, 5, 0.0100, 0.0850, 0.1760],
        [5, 7, 0.0320, 0.1610, 0.3060], [7, 2, 0.0000, 0.0625, 0.0000],
        [7, 8, 0.0085, 0.0720, 0.1490], [8, 9, 0.0119, 0.1008, 0.2090],
        [9, 3, 0.0000, 0.0586, 0.0000], [9, 6, 0.0390, 0.1700, 0.3580],
        [6, 4, 0.0170, 0.0920, 0.1580]
    ]
    line_data = [{'from': l[0], 'to': l[1], 'r': l[2], 'x': l[3], 'b': l[4]} for l in line_configs]
    
    return bus_data, line_data

def build_y_bus(bus_data, line_data):
    num_buses = len(bus_data)
    id_map = {b['id']: i for i, b in enumerate(bus_data)}
    Y = np.zeros((num_buses, num_buses), dtype=complex)
    for line in line_data:
        if line['from'] in id_map and line['to'] in id_map:
            i, j = id_map[line['from']], id_map[line['to']]
            z = complex(line['r'], line['x'])
            y_s = 1/z
            y_sh = complex(0, line['b']/2)
            Y[i, j] -= y_s
            Y[j, i] -= y_s
            Y[i, i] += (y_s + y_sh)
            Y[j, j] += (y_s + y_sh)
    return Y

# =================================================================
# 2. NEWTON-RAPHSON SOLVER
# =================================================================
def run_nr_loadflow(Y_bus, bus_data, max_iter=15, tol=1e-6):
    num_buses = len(bus_data)
    V = np.array([b['V'] for b in bus_data])
    Theta = np.array([b['theta'] for b in bus_data])
    P_spec = np.array([b['P_spec'] for b in bus_data])
    Q_spec = np.array([b['Q_spec'] for b in bus_data])
    types = np.array([b['type'] for b in bus_data]) 

    idx_slack = np.where(types == 1)[0]
    idx_pv    = np.where(types == 2)[0]
    idx_pq    = np.where(types == 3)[0]
    non_slack = np.sort(np.concatenate((idx_pv, idx_pq)))
    
    for it in range(max_iter):
        P_calc = np.zeros(num_buses)
        Q_calc = np.zeros(num_buses)
        for i in range(num_buses):
            for k in range(num_buses):
                mag_Y, ang_Y = abs(Y_bus[i, k]), np.angle(Y_bus[i, k])
                P_calc[i] += V[i] * V[k] * mag_Y * np.cos(Theta[i] - Theta[k] - ang_Y)
                Q_calc[i] += V[i] * V[k] * mag_Y * np.sin(Theta[i] - Theta[k] - ang_Y)

        dPa, dQa = P_spec - P_calc, Q_spec - Q_calc
        
        # Mismatch check
        mismatch = np.max(np.abs(dPa[non_slack]))
        if len(idx_pq) > 0:
            mismatch = max(mismatch, np.max(np.abs(dQa[idx_pq])))
            
        if mismatch < tol:
            return V, Theta, P_calc, Q_calc

        J1, J2, J3, J4 = [np.zeros((num_buses, num_buses)) for _ in range(4)]
        for i in range(num_buses):
            for k in range(num_buses):
                mag_Y, ang_Y = abs(Y_bus[i, k]), np.angle(Y_bus[i, k])
                if i != k:
                    J1[i, k] = V[i] * V[k] * mag_Y * np.sin(Theta[i] - Theta[k] - ang_Y)
                    J2[i, k] = V[i] * mag_Y * np.cos(Theta[i] - Theta[k] - ang_Y)
                    J3[i, k] = -V[i] * V[k] * mag_Y * np.cos(Theta[i] - Theta[k] - ang_Y)
                    J4[i, k] = V[i] * mag_Y * np.sin(Theta[i] - Theta[k] - ang_Y)
                else:
                    J1[i, i] = -Q_calc[i] - (V[i]**2 * Y_bus[i,i].imag)
                    J2[i, i] = P_calc[i] / V[i] + (V[i] * Y_bus[i,i].real)
                    J3[i, i] = P_calc[i] - (V[i]**2 * Y_bus[i,i].real)
                    J4[i, i] = Q_calc[i] / V[i] - (V[i] * Y_bus[i,i].imag)

        J = np.block([[J1[np.ix_(non_slack, non_slack)], J2[np.ix_(non_slack, idx_pq)]],
                      [J3[np.ix_(idx_pq, non_slack)], J4[np.ix_(idx_pq, idx_pq)]]])
        M = np.concatenate((dPa[non_slack], dQa[idx_pq]))
        
        try:
            corr = np.linalg.solve(J, M)
            Theta[non_slack] += corr[:len(non_slack)]
            V[idx_pq] += corr[len(non_slack):]
        except: return None, None, None, None

    return V, Theta, P_calc, Q_calc

# =================================================================
# 3. MAIN SIMULATION ENGINE
# =================================================================
def main():
    global SYSTEM_FREQ
    b_data, l_data = get_ieee_9bus_data()
    
    # 1. User selects generator to trip
    pv_buses = [b for b in b_data if b['type'] == 2]
    print("\n" + "="*40)
    print("      SELECT GENERATOR TO TRIP")
    print("="*40)
    for b in pv_buses:
        print(f"ID: {b['id']} | Pg: {b['Pg']:.4f} pu")
    target_trip_id = int(input("Enter Bus ID to trip: "))

    Y_bus = build_y_bus(b_data, l_data)
    
    # 2. Setup Real-time Plot
    plt.ion()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    times, freqs, rocofs = [], [], []
    
    ax1.set_title('Real-Time Frequency Response (IEEE 9-Bus)')
    ax1.set_ylabel('Frequency (Hz)')
    ax1.set_xlim(0, TOTAL_TIME)
    ax1.set_ylim(48.5, 50.5)
    line_f, = ax1.plot([], [], 'r-', label='System Frequency')
    ax1.grid(True)
    
    ax2.set_ylabel('RoCoF (Hz/s)')
    ax2.set_xlabel('Time (s)')
    ax2.set_xlim(0, TOTAL_TIME)
    ax2.set_ylim(-0.5, 0.1)
    line_r, = ax2.plot([], [], 'b-', label='RoCoF')
    ax2.grid(True)

    print(f"\n--- Starting Simulation (Trip at t={TRIP_TIME}s) ---")
    
    for t in np.arange(0, TOTAL_TIME + TIME_STEP, TIME_STEP):
        # --- Contingency Logic [cite: 14, 25] ---
        if abs(t - TRIP_TIME) < 1e-3:
            print(f"\n{RED}!!! EVENT: BUS {target_trip_id} TRIPPED !!!{RESET}")
            b_data = [b for b in b_data if b['id'] != target_trip_id]
            Y_bus = build_y_bus(b_data, l_data)

        # --- Solve Load Flow [cite: 23] ---
        V_sol, Th_sol, P_calc, Q_calc = run_nr_loadflow(Y_bus, b_data)
        
        if V_sol is None:
            print(f"{RED}Simulation Crash: Voltage Collapse at t={t}s{RESET}")
            break

        # --- Frequency Physics Logic ---
        # 1. Find Slack Generator Power (The generator that balances the system)
        net_imbalance = 0.0
        for i, b in enumerate(b_data):
            if b['type'] == 1: # Slack Bus
                p_required = P_calc[i] + b['Pl'] # P_calc is net injection (Pg - Pl)
                if p_required > b['P_max']:
                    # Deficit: What we need minus what we can physically provide
                    net_imbalance = b['P_max'] - p_required 
                
        # 2. Calculate Swing Equation [cite: 25]
        total_load = sum([b['Pl'] for b in b_data])
        rocof = 0.0
        if abs(net_imbalance) > 1e-6:
            # df/dt = (DeltaP * f_nom) / (2 * H * P_load)
            rocof = (net_imbalance * 50.0) / (2 * H_CONST * total_load)
            SYSTEM_FREQ += rocof * TIME_STEP
        
        # Update Plotting Data
        times.append(t)
        freqs.append(SYSTEM_FREQ)
        rocofs.append(rocof)
        
        line_f.set_data(times, freqs)
        line_r.set_data(times, rocofs)
        plt.pause(0.05)
        
        # Update bus states for next iteration
        for i in range(len(b_data)):
            b_data[i]['V'] = V_sol[i]
            b_data[i]['theta'] = Th_sol[i]

        print(f"t={t:4.1f}s | Freq: {SYSTEM_FREQ:.4f} Hz | RoCoF: {rocof:7.4f} Hz/s")

    plt.ioff()
    print("\nSimulation Complete.")
    plt.show()

if __name__ == "__main__":
    main()
