import tkinter as tk
from tkinter import scrolledtext, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Constants
G = 6.67430e-11  # Gravitational constant
M = 5.972e24  # Mass of Earth in kg
R_EARTH = 6371  # Radius of Earth in km

def parse_tle(tle_line1, tle_line2):
    """Parses TLE lines to extract orbital elements."""
    try:
        # Extracting from line 1
        inclination = float(tle_line2[8:16])
        raan = float(tle_line2[17:25])
        eccentricity = float("." + tle_line2[26:33])
        arg_perigee = float(tle_line2[34:42])
        mean_anomaly = float(tle_line2[43:51])
        mean_motion = float(tle_line2[52:63])
        return {
            "inclination": inclination,
            "raan": raan,
            "eccentricity": eccentricity,
            "arg_perigee": arg_perigee,
            "mean_anomaly": mean_anomaly,
            "mean_motion": mean_motion,
        }
    except (ValueError, IndexError):
        return None

def calculate_orbital_metrics(elements):
    """Calculates derived orbital metrics from orbital elements."""
    if not elements:
        return None
    try:
        mean_motion_rad_s = elements["mean_motion"] * 2 * np.pi / 86400
        semi_major_axis_m = ((G * M) / mean_motion_rad_s**2)**(1/3)
        semi_major_axis_km = semi_major_axis_m / 1000
        apogee = semi_major_axis_km * (1 + elements["eccentricity"]) - R_EARTH
        perigee = semi_major_axis_km * (1 - elements["eccentricity"]) - R_EARTH
        orbital_period = (2 * np.pi) / mean_motion_rad_s / 60  # in minutes
        return {
            "semi_major_axis": semi_major_axis_km,
            "apogee": apogee,
            "perigee": perigee,
            "orbital_period": orbital_period,
        }
    except (ValueError, ZeroDivisionError):
        return None

def estimate_visibility(elements, lat, lon):
    """Provides a simple visibility estimate."""
    if not elements or lat is None or lon is None:
        return "Observer location or satellite data missing."
    # Simplified logic: check if inclination is greater than observer's latitude
    if elements["inclination"] > abs(lat):
        return "Satellite may be visible."
    else:
        return "Satellite unlikely to be visible."

def plot_orbit(elements):
    """Plots a 2D projection of the satellite orbit."""
    fig, ax = plt.subplots(figsize=(5, 5))
    if elements:
        a = calculate_orbital_metrics(elements)["semi_major_axis"]
        e = elements["eccentricity"]
        b = a * np.sqrt(1 - e**2)
        
        # Parametric equation for an ellipse
        t = np.linspace(0, 2 * np.pi, 100)
        x = a * np.cos(t)
        y = b * np.sin(t)
        
        # Rotation for argument of perigee
        arg_p_rad = np.deg2rad(elements["arg_perigee"])
        x_rot = x * np.cos(arg_p_rad) - y * np.sin(arg_p_rad)
        y_rot = x * np.sin(arg_p_rad) + y * np.cos(arg_p_rad)
        
        ax.plot(x_rot, y_rot, label="Satellite Orbit")
    
    # Plot Earth
    earth = plt.Circle((0, 0), R_EARTH, color='blue', label='Earth')
    ax.add_artist(earth)
    
    ax.set_title("2D Satellite Orbit Projection")
    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True)
    
    # Set limits to show the full orbit
    if elements:
        max_dim = calculate_orbital_metrics(elements)["semi_major_axis"] * 1.2
        ax.set_xlim(-max_dim, max_dim)
        ax.set_ylim(-max_dim, max_dim)
    else:
        ax.set_xlim(-10000, 10000)
        ax.set_ylim(-10000, 10000)
        
    return fig

class SatelliteApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Satellite Service UI")
        
        # --- Input Fields ---
        tk.Label(root, text="TLE Line 1:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.tle1_entry = tk.Entry(root, width=80)
        self.tle1_entry.grid(row=0, column=1, padx=5, pady=5)
        
        tk.Label(root, text="TLE Line 2:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.tle2_entry = tk.Entry(root, width=80)
        self.tle2_entry.grid(row=1, column=1, padx=5, pady=5)
        
        tk.Label(root, text="Observer Latitude:").grid(row=2, column=0, padx=5, pady=5, sticky="w")
        self.lat_entry = tk.Entry(root, width=20)
        self.lat_entry.grid(row=2, column=1, padx=5, pady=5, sticky="w")
        
        tk.Label(root, text="Observer Longitude:").grid(row=3, column=0, padx=5, pady=5, sticky="w")
        self.lon_entry = tk.Entry(root, width=20)
        self.lon_entry.grid(row=3, column=1, padx=5, pady=5, sticky="w")
        
        # --- Buttons ---
        self.process_button = tk.Button(root, text="Process Data", command=self.process_data)
        self.process_button.grid(row=4, column=0, padx=5, pady=10)
        
        self.clear_button = tk.Button(root, text="Clear Inputs", command=self.clear_inputs)
        self.clear_button.grid(row=4, column=1, padx=5, pady=10, sticky="w")
        
        # --- Results Area ---
        self.results_text = scrolledtext.ScrolledText(root, width=100, height=10)
        self.results_text.grid(row=5, column=0, columnspan=2, padx=5, pady=5)
        
        # --- Plot Area ---
        self.fig = plot_orbit(None)
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().grid(row=6, column=0, columnspan=2, padx=5, pady=5)
        self.canvas.draw()

    def process_data(self):
        tle1 = self.tle1_entry.get()
        tle2 = self.tle2_entry.get()
        lat_str = self.lat_entry.get()
        lon_str = self.lon_entry.get()
        
        if not tle1 or not tle2 or not lat_str or not lon_str:
            messagebox.showerror("Input Error", "All fields are required.")
            return
        
        try:
            lat = float(lat_str)
            lon = float(lon_str)
        except ValueError:
            messagebox.showerror("Input Error", "Latitude and Longitude must be numbers.")
            return
            
        elements = parse_tle(tle1, tle2)
        if not elements:
            messagebox.showerror("Input Error", "Invalid TLE data.")
            return
            
        metrics = calculate_orbital_metrics(elements)
        if not metrics:
            messagebox.showerror("Processing Error", "Could not compute orbital metrics.")
            return
            
        visibility = estimate_visibility(elements, lat, lon)
        
        # Display results
        self.results_text.delete(1.0, tk.END)
        self.results_text.insert(tk.END, "--- Orbital Metrics ---\n")
        self.results_text.insert(tk.END, f"Semi-Major Axis: {metrics['semi_major_axis']:.2f} km\n")
        self.results_text.insert(tk.END, f"Apogee: {metrics['apogee']:.2f} km\n")
        self.results_text.insert(tk.END, f"Perigee: {metrics['perigee']:.2f} km\n")
        self.results_text.insert(tk.END, f"Orbital Period: {metrics['orbital_period']:.2f} minutes\n\n")
        self.results_text.insert(tk.END, "--- Visibility Estimate ---\n")
        self.results_text.insert(tk.END, f"{visibility}\n")
        
        # Update plot
        self.fig = plot_orbit(elements)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().grid(row=6, column=0, columnspan=2, padx=5, pady=5)
        self.canvas.draw()

    def clear_inputs(self):
        self.tle1_entry.delete(0, tk.END)
        self.tle2_entry.delete(0, tk.END)
        self.lat_entry.delete(0, tk.END)
        self.lon_entry.delete(0, tk.END)
        self.results_text.delete(1.0, tk.END)
        
        # Clear plot
        self.fig = plot_orbit(None)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().grid(row=6, column=0, columnspan=2, padx=5, pady=5)
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = SatelliteApp(root)
    root.mainloop()