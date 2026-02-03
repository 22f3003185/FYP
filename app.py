import numpy as np
import cmath
import copy
from flask import Flask, jsonify, render_template, request
from flask_cors import CORS

app = Flask(__name__, template_folder='.')
CORS(app)

# --- GEN DATA (Inertia H, Transient Reactance Xd') ---
GEN_DATA = {
    0: {'H': 23.64, 'Xd_p': 0.0608}, # Gen 1
    1: {'H': 6.40,  'Xd_p': 0.1198}, # Gen 2
    2: {'H': 3.01,  'Xd_p': 0.1813}  # Gen 3
}

class PowerSystemSolver:
    def __init__(self):
        self.n_buses = 9
        self.freq = 50.0
        self.omega_s = 2 * np.pi * self.freq
        
        # [Bus, Type, V_spec, Theta, P_gen, Q_gen, P_load, Q_load]
        self.bus_data = np.array([
            [1, 1, 1.040, 0.0, 0.00, 0.0, 0.00, 0.00],
            [2, 2, 1.025, 0.0, 1.63, 0.0, 0.00, 0.00],
            [3, 2, 1.025, 0.0, 0.85, 0.0, 0.00, 0.00],
            [4, 3, 1.000, 0.0, 0.00, 0.0, 0.00, 0.00],
            [5, 3, 1.000, 0.0, 0.00, 0.0, 1.25, 0.50],
            [6, 3, 1.000, 0.0, 0.00, 0.0, 0.90, 0.30],
            [7, 3, 1.000, 0.0, 0.00, 0.0, 0.00, 0.00],
            [8, 3, 1.000, 0.0, 0.00, 0.0, 1.00, 0.35],
            [9, 3, 1.000, 0.0, 0.00, 0.0, 0.00, 0.00]
        ])
        
        # [From, To, R, X, B]
        self.branch_data = np.array([
            [1, 4, 0.0000, 0.0576, 0.0000], [4, 5, 0.0100, 0.0850, 0.1760],
            [5, 7, 0.0320, 0.1610, 0.3060], [7, 2, 0.0000, 0.0625, 0.0000],
            [7, 8, 0.0085, 0.0720, 0.1490], [8, 9, 0.0119, 0.1008, 0.2090],
            [9, 3, 0.0000, 0.0586, 0.0000], [9, 6, 0.0390, 0.1700, 0.3580],
            [6, 4, 0.0170, 0.0920, 0.1580]
        ])

        self.V = self.bus_data[:, 2].copy()
        self.delta = np.zeros(self.n_buses)
        self.Y_bus = self.build_ybus(self.branch_data)

    def build_ybus(self, branches):
        Y = np.zeros((self.n_buses, self.n_buses), dtype=complex)
        for br in branches:
            f, t = int(br[0])-1, int(br[1])-1
            z = complex(br[2], br[3])
            y_s = 1/z if abs(z)>0 else 1000
            y_sh = complex(0, br[4]/2)
            Y[f, t] -= y_s
            Y[t, f] -= y_s
            Y[f, f] += y_s + y_sh
            Y[t, t] += y_s + y_sh
        return Y

    def solve_load_flow(self):
        # Reset to flat start for clean solution
        self.V = self.bus_data[:, 2].copy()
        self.delta = np.zeros(self.n_buses)
        
        # Standard Newton-Raphson
        for _ in range(15):
            P_calc = np.zeros(9); Q_calc = np.zeros(9)
            for i in range(9):
                for k in range(9):
                    ang = self.delta[i] - self.delta[k]
                    term = self.V[i]*self.V[k]*abs(self.Y_bus[i,k])
                    theta = cmath.phase(self.Y_bus[i,k])
                    P_calc[i] += term * np.cos(theta - ang)
                    Q_calc[i] -= term * np.sin(theta - ang)
            
            for i in range(1, 9):
                P_spec = self.bus_data[i, 4] - self.bus_data[i, 6]
                dP = P_spec - P_calc[i]
                # dP/dDelta approx -B * V^2
                J_p = -self.Y_bus[i,i].imag * self.V[i]**2
                if abs(J_p) > 0.0001: self.delta[i] += dP / J_p
                
                if self.bus_data[i, 1] == 3: # PQ Bus
                    Q_spec = self.bus_data[i, 5] - self.bus_data[i, 7]
                    dQ = Q_spec - Q_calc[i]
                    # dQ/dV approx -B * 2V
                    J_q = -self.Y_bus[i,i].imag * 2 * self.V[i]
                    if abs(J_q) > 0.0001: self.V[i] += dQ / J_q
                    
        self.P_net = P_calc
        self.Q_net = Q_calc

    def simulate_fault(self, fault_line_id):
        # 1. Start from Steady State
        self.solve_load_flow()
        
        # 2. Init Generator Models
        gen_state = {}
        for idx in GEN_DATA:
            S = complex(self.P_net[idx], self.Q_net[idx])
            V_p = self.V[idx] * cmath.exp(1j * self.delta[idx])
            I_inj = np.conj(S / V_p)
            E_prime = V_p + (1j * GEN_DATA[idx]['Xd_p'] * I_inj)
            gen_state[idx] = {'E': abs(E_prime), 'delta': cmath.phase(E_prime), 'w': self.omega_s, 'Pm': self.P_net[idx]}

        # 3. Build Fault Y-Bus
        # Map frontend ID to bus index (simple mapping)
        map_ids = {"line_4_5": 4, "line_4_6": 4, "line_5_7": 5, "line_6_9": 6, "line_7_8": 7, "line_8_9": 8}
        fault_bus = map_ids.get(fault_line_id, 4) - 1
        
        Y_f = copy.deepcopy(self.Y_bus)
        Y_f[fault_bus, fault_bus] += complex(0, -1e6) # Short Circuit
        
        # Add Load & Gen Admittances
        for i in range(9):
            if self.bus_data[i, 1] == 3:
                S_load = complex(self.bus_data[i, 6], self.bus_data[i, 7])
                Y_f[i, i] += np.conj(S_load) / (self.V[i]**2)
        
        for idx in GEN_DATA:
            Y_f[idx, idx] += 1 / (1j * GEN_DATA[idx]['Xd_p'])
            
        Z_f = np.linalg.inv(Y_f)

        # 4. Integrate (0.1s)
        dt = 0.005
        for _ in range(int(0.1/dt)):
            # I_source injection
            I_vec = np.zeros(9, dtype=complex)
            for idx in GEN_DATA:
                E = gen_state[idx]['E'] * cmath.exp(1j * gen_state[idx]['delta'])
                I_vec[idx] = E / (1j * GEN_DATA[idx]['Xd_p'])
            
            V_cplx = np.dot(Z_f, I_vec)
            
            # Update Gen State
            for idx in GEN_DATA:
                E = gen_state[idx]['E'] * cmath.exp(1j * gen_state[idx]['delta'])
                Pe = (E * np.conj((E - V_cplx[idx])/(1j * GEN_DATA[idx]['Xd_p']))).real
                
                dw = (np.pi * 50 / GEN_DATA[idx]['H']) * (gen_state[idx]['Pm'] - Pe)
                gen_state[idx]['w'] += dw * dt
                gen_state[idx]['delta'] += (gen_state[idx]['w'] - self.omega_s) * dt
        
        # Save Fault State
        self.V = np.abs(V_cplx)
        self.delta = np.angle(V_cplx)
        # Recalc Powers for display
        P, Q = np.zeros(9), np.zeros(9) # Approximate for display
        self.P_net = P; self.Q_net = Q 

    def get_data(self):
        res = {"buses": {}, "lines": {}}
        
        # Buses
        for i in range(9):
            bid = f"bus_{i+1}"
            res["buses"][bid] = {
                "name": f"Bus {i+1}",
                "pu": round(self.V[i], 4),
                "kv": round(self.V[i] * (16.5 if i==0 else 18 if i==1 else 13.8 if i==2 else 230), 2),
                "angle": round(np.degrees(self.delta[i]), 2),
                "p": round(self.P_net[i]*100, 2),
                "q": round(self.Q_net[i]*100, 2),
                "freq": 50.0
            }
            
        # Lines
        topo = {
            (1,4):"tf_4_1", (4,1):"tf_4_1", (4,5):"line_4_5", (5,4):"line_4_5",
            (4,6):"line_4_6", (6,4):"line_4_6", (5,7):"line_5_7", (7,5):"line_5_7",
            (6,9):"line_6_9", (9,6):"line_6_9", (7,8):"line_7_8", (8,7):"line_7_8",
            (8,9):"line_8_9", (9,8):"line_8_9", (2,7):"tf_7_2", (7,2):"tf_7_2",
            (3,9):"tf_9_3", (9,3):"tf_9_3"
        }
        
        for br in self.branch_data:
            f, t = int(br[0]), int(br[1])
            fid = topo.get((f,t))
            if fid:
                vf = self.V[f-1]*cmath.exp(1j*self.delta[f-1])
                vt = self.V[t-1]*cmath.exp(1j*self.delta[t-1])
                y = 1/complex(br[2], br[3])
                s = vf * np.conj((vf-vt)*y) * 100
                res["lines"][fid] = {
                    "name": fid.upper(), "from":f"Bus {f}", "to":f"Bus {t}",
                    "mw": round(s.real, 2), "mvar": round(s.imag, 2), "mva": round(abs(s), 2)
                }
        return res

solver = PowerSystemSolver()

@app.route('/')
def home(): return render_template('index.html')

@app.route('/api/solve', methods=['POST'])
def solve():
    data = request.json
    if data.get('type') == 'fault':
        solver.simulate_fault(data.get('line_id'))
    else:
        solver.solve_load_flow()
    return jsonify(solver.get_data())

if __name__ == '__main__':
    app.run(debug=True, port=5000)