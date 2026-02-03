import numpy as np
import cmath
import copy
from flask import Flask, jsonify, render_template, request
from flask_cors import CORS

app = Flask(__name__, template_folder='.')
CORS(app)

# --- GEN DATA (Inertia H, Transient Reactance Xd') ---
GEN_DATA = {
    0: {'H': 23.64, 'Xd_p': 0.0608},
    1: {'H': 6.40,  'Xd_p': 0.1198},
    2: {'H': 3.01,  'Xd_p': 0.1813}
}

BASE_FREQ = 50.0
OMEGA_S = 2 * np.pi * BASE_FREQ


class PowerSystemSolver:
    def __init__(self):
        self.n_buses = 9

        # [Bus, Type, V, θ, Pg, Qg, Pl, Ql]
        self.bus_data = np.array([
            [1,1,1.040,0,0,0,0,0],
            [2,2,1.025,0,1.63,0,0,0],
            [3,2,1.025,0,0.85,0,0,0],
            [4,3,1.000,0,0,0,0,0],
            [5,3,1.000,0,0,0,1.25,0.50],
            [6,3,1.000,0,0,0,0.90,0.30],
            [7,3,1.000,0,0,0,0,0],
            [8,3,1.000,0,0,0,1.00,0.35],
            [9,3,1.000,0,0,0,0,0]
        ])

        self.branch_data = np.array([
            [1,4,0,0.0576,0],[4,5,0.01,0.085,0.176],
            [5,7,0.032,0.161,0.306],[7,2,0,0.0625,0],
            [7,8,0.0085,0.072,0.149],[8,9,0.0119,0.1008,0.209],
            [9,3,0,0.0586,0],[9,6,0.039,0.17,0.358],
            [6,4,0.017,0.092,0.158]
        ])

        self.Y_bus = self.build_ybus(self.branch_data)
        self.reset_steady_state()

        # ---- DYNAMIC STATE ----
        self.gen_state = {}
        self.stage = "steady"
        self.history = {"time": [], "freq": {}, "rocof": {}}

    # ------------------------------------------------------------------

    def reset_steady_state(self):
        self.V = self.bus_data[:,2].copy()
        self.delta = np.zeros(self.n_buses)

    def build_ybus(self, branches):
        Y = np.zeros((9,9), dtype=complex)
        for br in branches:
            f,t = int(br[0])-1, int(br[1])-1
            z = complex(br[2], br[3])
            y = 1/z if abs(z)>0 else 1e6
            ysh = complex(0, br[4]/2)
            Y[f,f]+=y+ysh; Y[t,t]+=y+ysh
            Y[f,t]-=y; Y[t,f]-=y
        return Y

    # ------------------------------------------------------------------

    def solve_load_flow(self):
        self.reset_steady_state()
        for _ in range(15):
            P,Q = np.zeros(9),np.zeros(9)
            for i in range(9):
                for k in range(9):
                    ang = self.delta[i]-self.delta[k]
                    Yik = self.Y_bus[i,k]
                    P[i]+=self.V[i]*self.V[k]*abs(Yik)*np.cos(cmath.phase(Yik)-ang)
                    Q[i]-=self.V[i]*self.V[k]*abs(Yik)*np.sin(cmath.phase(Yik)-ang)

            for i in range(1,9):
                dP=(self.bus_data[i,4]-self.bus_data[i,6])-P[i]
                self.delta[i]+=dP/(-self.Y_bus[i,i].imag*self.V[i]**2)
                if self.bus_data[i,1]==3:
                    dQ=(self.bus_data[i,5]-self.bus_data[i,7])-Q[i]
                    self.V[i]+=dQ/(-2*self.Y_bus[i,i].imag*self.V[i])

        self.P_net,self.Q_net=P,Q
        self.init_generators()

    # ------------------------------------------------------------------

    def init_generators(self):
        self.gen_state={}
        for idx in GEN_DATA:
            Vp=self.V[idx]*cmath.exp(1j*self.delta[idx])
            S=complex(self.P_net[idx],self.Q_net[idx])
            I=np.conj(S/Vp)
            E=Vp+1j*GEN_DATA[idx]['Xd_p']*I
            self.gen_state[idx]={
                "E":abs(E),
                "delta":cmath.phase(E),
                "omega":OMEGA_S,
                "Pm":self.P_net[idx],
                "online":True
            }

    # ------------------------------------------------------------------

    def integrate(self,Y,dt,steps):
        self.history={"time":[],"freq":{},"rocof":{}}
        for g in self.gen_state:
            self.history["freq"][g]=[]
            self.history["rocof"][g]=[]

        Z=np.linalg.inv(Y)

        for k in range(steps):
            self.history["time"].append(k*dt)
            I=np.zeros(9,dtype=complex)

            for g in self.gen_state:
                if self.gen_state[g]["online"]:
                    E=self.gen_state[g]["E"]*cmath.exp(1j*self.gen_state[g]["delta"])
                    I[g]=E/(1j*GEN_DATA[g]["Xd_p"])

            V=np.dot(Z,I)

            for g in self.gen_state:
                if not self.gen_state[g]["online"]:
                    continue
                E=self.gen_state[g]["E"]*cmath.exp(1j*self.gen_state[g]["delta"])
                Pe=(E*np.conj((E-V[g])/(1j*GEN_DATA[g]["Xd_p"])) ).real
                H=GEN_DATA[g]["H"]
                domega=(np.pi*BASE_FREQ/H)*(self.gen_state[g]["Pm"]-Pe)
                self.gen_state[g]["omega"]+=domega*dt
                self.gen_state[g]["delta"]+=(self.gen_state[g]["omega"]-OMEGA_S)*dt

                f=self.gen_state[g]["omega"]/(2*np.pi)
                self.history["freq"][g].append(f)
                self.history["rocof"][g].append(domega/(2*np.pi))

        self.V=np.abs(V); self.delta=np.angle(V)

    # ------------------------------------------------------------------

    def simulate_fault(self,line_id):
        self.solve_load_flow()
        Yf=copy.deepcopy(self.Y_bus)
        bus_map={"line_4_5":4,"line_4_6":4,"line_5_7":5,"line_6_9":6,"line_7_8":7,"line_8_9":8}
        fb=bus_map.get(line_id,4)-1
        Yf[fb,fb]+=complex(0,-1e6)

        for i in range(9):
            if self.bus_data[i,1]==3:
                S=complex(self.bus_data[i,6],self.bus_data[i,7])
                Yf[i,i]+=np.conj(S)/(self.V[i]**2)

        for g in GEN_DATA:
            Yf[g,g]+=1/(1j*GEN_DATA[g]['Xd_p'])

        self.stage="fault"
        self.integrate(Yf,0.005,int(0.1/0.005))

    # ------------------------------------------------------------------

    def clear_fault(self):
        self.stage="postfault"
        self.integrate(self.Y_bus,0.005,int(0.2/0.005))

    # ------------------------------------------------------------------

    def trip_generator(self,gid):
        if gid in self.gen_state:
            self.gen_state[gid]["online"]=False
            self.gen_state[gid]["Pm"]=0.0
            self.clear_fault()

    # ------------------------------------------------------------------

    def shed_load(self,bus,scale):
        self.bus_data[bus,6]*=scale
        self.bus_data[bus,7]*=scale
        self.solve_load_flow()

    # ------------------------------------------------------------------

    def get_data(self):
        res = {"buses": {}, "lines": {}, "history": self.history}
        # Compute per-bus frequency and ROCOF (mean of online generators)
        gen_freqs = {}
        gen_rocofs = {}
        for g in self.gen_state:
            if self.gen_state[g]["online"]:
                gen_freqs[g] = self.gen_state[g]["omega"] / (2 * np.pi)
                # ROCOF: use last value from history if available, else 0
                if self.history.get("rocof", {}).get(g):
                    gen_rocofs[g] = self.history["rocof"][g][-1]
                else:
                    gen_rocofs[g] = 0.0
        mean_freq = np.mean(list(gen_freqs.values())) if gen_freqs else BASE_FREQ
        mean_rocof = np.mean(list(gen_rocofs.values())) if gen_rocofs else 0.0

        # Add bus fields: name, kv, p, q, rocof
        # For this system, assume all buses at 230kV, except buses 7-9 at 13.8kV (standard IEEE 9-bus)
        bus_kv = [16.5, 18.0, 13.8, 230.0, 230.0, 230.0, 230.0, 230.0, 230.0]
        for i in range(9):
            # Approximate bus frequency and rocof as mean of generator values
            bus_freq = mean_freq
            bus_rocof = mean_rocof
            # Bus name
            bus_name = f"Bus {i+1}"
            # Bus kV
            kv = bus_kv[i]
            # p, q (net injection, positive for generation, negative for load)
            p = round(self.P_net[i], 4) if hasattr(self, "P_net") else 0.0
            q = round(self.Q_net[i], 4) if hasattr(self, "Q_net") else 0.0
            res["buses"][f"bus_{i+1}"] = {
                "pu": round(self.V[i], 4),
                "angle": round(np.degrees(self.delta[i]), 2),
                "freq": round(bus_freq, 3),
                "name": bus_name,
                "kv": kv,
                "p": p,
                "q": q,
                "rocof": round(bus_rocof, 6),
            }

        # Compute line flows using branch_data and Y_bus
        # Use bus voltages (V, delta) to compute S_ij for each branch
        # Map line indices to SVG IDs
        branch_svg = [
            "line_1_4", "line_4_5", "line_5_7", "tf_7_2", "line_7_8",
            "line_8_9", "tf_9_3", "line_9_6", "line_6_4"
        ]
        for idx, br in enumerate(self.branch_data):
            f, t = int(br[0]) - 1, int(br[1]) - 1
            r, x, b = br[2], br[3], br[4]
            z = complex(r, x)
            y = 1 / z if abs(z) > 0 else 1e6
            ysh = complex(0, b)
            V_f = self.V[f] * cmath.exp(1j * self.delta[f])
            V_t = self.V[t] * cmath.exp(1j * self.delta[t])
            # Line current from f to t
            I_ft = (V_f - V_t) * y + V_f * ysh / 2
            S_ft = V_f * np.conj(I_ft)
            # Line current from t to f
            I_tf = (V_t - V_f) * y + V_t * ysh / 2
            S_tf = V_t * np.conj(I_tf)
            # SVG ID: handle transformers (tf_7_2, tf_9_3)
            svg_id = branch_svg[idx]
            res["lines"][svg_id] = {
                "from": f"bus_{f+1}",
                "to": f"bus_{t+1}",
                "p_from": round(S_ft.real, 4),
                "q_from": round(S_ft.imag, 4),
                "p_to": round(S_tf.real, 4),
                "q_to": round(S_tf.imag, 4),
            }

        return res


solver=PowerSystemSolver()

@app.route("/")
def home(): return render_template("index_def.html")

@app.route("/api/solve",methods=["POST"])
def solve():
    d=request.json
    t=d.get("type")
    if t=="fault": solver.simulate_fault(d.get("line_id"))
    elif t=="postfault": solver.clear_fault()
    elif t=="gen_trip": solver.trip_generator(int(d.get("gen_id")[-1])-1)
    elif t=="load_adjust": solver.shed_load(int(d.get("load_id")[-1])-1,d["scale"])
    else: solver.solve_load_flow()
    return jsonify(solver.get_data())

if __name__=="__main__":
    app.run(debug=True,port=5000)