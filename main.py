import numpy as np
import plotly.graph_objects as go
from nicegui import app, ui



# === 1. DEFINE YOUR SIMULATION FUNCTION (CORRECTED) ===
def calculate_trajectory(g, m, theta_deg, v0, c_d, r, rho_air, dt, air_resistance_on: bool):
    """
    Runs the projectile motion simulation and returns a Plotly figure.
    'g' is now a variable.
    """
    
    # --- Define Constants ---
    A = np.pi * r**2  # Cross-sectional area (m^2)

    # --- Initialize Lists ---
    x_points = []
    y_points = []
    
    # --- Initialize Plotly Figure ---
    fig = go.Figure()
    
    theta_rad = np.radians(theta_deg)

    # --- BRANCH 1: CALCULATE WITH AIR RESISTANCE ---
    if air_resistance_on:
        # --- Initialize Simulation ---
        x = 0.0
        y = 0.0
        v_x = v0 * np.cos(theta_rad)
        v_y = v0 * np.sin(theta_rad)
        x_points = [x]
        y_points = [y]

        # --- Run the Simulation Loop (The Drag Calculation) ---
        while y >= 0:
            v = np.sqrt(v_x**2 + v_y**2)
            
            if v == 0:
                F_drag_x, F_drag_y = 0, 0
            else:
                F_air = 0.5 * c_d * rho_air * A * (v**2)
                F_drag_x = -F_air * (v_x / v)
                F_drag_y = -F_air * (v_y / v)
            
            F_g = -m * g 
            F_total_x = F_drag_x
            F_total_y = F_g + F_drag_y

            a_x = F_total_x / m
            a_y = F_total_y / m
            
            v_x += (a_x * dt)
            v_y += (a_y * dt)
            x += (v_x * dt)
            y += (v_y * dt)
            
            if y >= 0:
                x_points.append(x)
                y_points.append(y)
                
        fig.add_trace(go.Scatter(
            x=x_points, 
            y=y_points, 
            mode='lines', 
            name='Gravity and Drag',
            line=dict(color='blue', width=3)
        ))

    # --- BRANCH 2: CALCULATE IDEAL PATH (Gravity Only) ---
    v_x0_ideal = v0 * np.cos(theta_rad)
    v_y0_ideal = v0 * np.sin(theta_rad)
    
    if g == 0:
        t_f_ideal = 0 
    elif v_y0_ideal < 0:
        t_f_ideal = 0 
    else:
        t_f_ideal = (2 * v_y0_ideal) / g

    num_points = len(x_points) if air_resistance_on and len(x_points) > 0 else 200
    
    if t_f_ideal == 0:
         t_ideal = np.zeros(num_points)
    else:
         t_ideal = np.linspace(0, t_f_ideal, num_points)

    x_ideal = v_x0_ideal * t_ideal
    y_ideal = (v_y0_ideal * t_ideal) - (0.5 * g * t_ideal**2)

    fig.add_trace(go.Scatter(
        x=x_ideal, 
        y=y_ideal, 
        mode='lines', 
        name='Gravity Only (Ideal)',
        line=dict(color='black', dash='dot')
    ))

    # --- Update layout ---
    title_suffix = "With vs. Without Air Resistance" if air_resistance_on else "Gravity Only (Ideal)"
    
    # --- THIS IS THE FIX ---
    # Safely get max values from lists (if they exist)
    y_max_drag = max(y_points) if y_points else 0
    x_max_drag = max(x_points) if x_points else 0
    
    # Safely get max values from numpy arrays (if they exist)
    y_max_ideal = y_ideal.max() if y_ideal.size > 0 else 0
    x_max_ideal = x_ideal.max() if x_ideal.size > 0 else 0

    # Now find the overall max
    yaxis_max = max(y_max_drag, y_max_ideal) * 1.1 + 1
    xaxis_max = max(x_max_drag, x_max_ideal) * 1.1 + 1
    # --- END OF FIX ---

    fig.update_layout(
        title=f'Trajectory (v0={v0} m/s, angle={theta_deg}°, g={g} m/s²) - {title_suffix}',
        xaxis_title='Displacement in x-direction (m)',
        yaxis_title='Displacement in y-direction (m)',
        legend_title='Simulation Model',
        
        yaxis=dict(scaleanchor="x", scaleratio=1, range=[0, yaxis_max]),
        xaxis=dict(range=[0, xaxis_max]),
        margin=dict(l=20, r=20, t=40, b=20)
    )
    
    return fig


# === 2. CREATE THE NICEGUI UI ===
# (This is your existing main_page function - it's correct!)
@ui.page('/')
def main_page():
    
    # --- 1. Define Header ---
    with ui.header(elevated=True).classes('items-center justify-between px-4 bg-gray-800 text-white') as header:
        with ui.row().classes() as menu_row:
             ui.button(on_click=lambda: left_drawer.toggle(), icon='menu').props('flat color=white').classes('mr-2')
             ui.label('Projectile Motion Simulator').classes('text-lg font-bold')
             

    # --- 2. Define Left Drawer ---
    with ui.left_drawer(elevated=True).classes('bg-gray-100 p-4') as left_drawer:
        ui.label('Simulation Parameters').classes('text-lg font-bold mb-2')
        ui.separator()

        # --- We now store the inputs in a dictionary ---
        inputs = {}
        inputs['g'] = ui.number(label='Gravity (m/s²)', value=9.8, min=0, step=0.1, format='%.1f').props('outlined dense')
        inputs['m'] = ui.number(label='Mass (kg)', value=2.0, min=0.1, step=0.1, format='%.1f').props('outlined dense')
        inputs['theta_deg'] = ui.number(label='Angle (degrees)', value=45.0, min=0, max=90, step=1, format='%.0f').props('outlined dense')
        inputs['v0'] = ui.number(label='Initial Velocity (m/s)', value=100.0, min=1, step=1, format='%.0f').props('outlined dense')
        
        inputs['air_resistance'] = ui.switch('Air Resistance', value=True).classes('mt-2')
        
        inputs['c_d'] = ui.number(label='Drag Coefficient', value=0.45, min=0, step=0.01, format='%.2f').props('outlined dense')
        inputs['r'] = ui.number(label='Radius (m)', value=0.5, min=0.01, step=0.01, format='%.2f').props('outlined dense')
        inputs['rho_air'] = ui.number(label='Air Density (kg/m^3)', value=0.0175, min=0.001, step=0.001, format='%.4f').props('outlined dense')
        
        # Bind the inputs to the switch
        inputs['c_d'].bind_enabled_from(inputs['air_resistance'], 'value')
        inputs['r'].bind_enabled_from(inputs['air_resistance'], 'value')
        inputs['rho_air'].bind_enabled_from(inputs['air_resistance'], 'value')

        inputs['dt'] = ui.number(label='Time Step (s)', value=0.01, min=0.001, step=0.001, format='%.3f').props('outlined dense')

        ui.separator().classes('my-4')
        
        # --- Create a placeholder for the button ---
        button_placeholder = ui.row().classes('w-full')


    # --- 4. Define Main Content Area ---
    with ui.column().classes('w-full h-screen p-4 items-center justify-center'):
        with ui.card().classes('w-full max-w-5xl h-[80vh]'):
            # --- Store the plot_view in a variable ---
            plot_view = ui.plotly(go.Figure())

    # --- 5. Define the update_plot function ---
    # It is defined *after* 'inputs' and 'plot_view' exist.
    def update_plot():
        try:
            # Read all values from the 'inputs' dictionary
            g = inputs['g'].value
            m = inputs['m'].value
            theta_deg = inputs['theta_deg'].value
            v0 = inputs['v0'].value
            c_d = inputs['c_d'].value
            r = inputs['r'].value
            rho_air = inputs['rho_air'].value
            dt = inputs['dt'].value
            air_resistance_on = inputs['air_resistance'].value

            # Call the calculation function
            fig = calculate_trajectory(g, m, theta_deg, v0, c_d, r, rho_air, dt, air_resistance_on)
            
            # Update the plot
            plot_view.update_figure(fig)
            
            ui.notify('Simulation complete!', type='positive')
        except Exception as e:
            ui.notify(f'Error during calculation: {e}', type='negative')
            print(e) 

    # --- 6. Add the button ---
    # Now that 'update_plot' exists, we can create the button
    with button_placeholder:
        ui.button('Run Simulation', on_click=update_plot, icon='play_arrow').classes('w-full')

    # --- 7. Initial Plot Load ---
    # Call the function once at the start
    update_plot()


# === 3. Run the app ===
ui.run(title="Projectile motion",favicon="🚀",port=8080,reload=False,storage_secret="superkey",host="0.0.0.0")
