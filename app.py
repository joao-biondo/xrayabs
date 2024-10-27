from flask import Flask, render_template, request, jsonify, flash, redirect
import xraydb as xr
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import re
import scipy.constants
import math, plotly

app = Flask(__name__)
app.secret_key = 'my_secret_key'

h = scipy.constants.physical_constants['Planck constant in eV/Hz'][0]
c = scipy.constants.c
mu = scipy.constants.m_u * (1e3)

def get_elements(chemical_formula):
    pattern = r'([A-Z][a-z]?)(\d*(\.\d+)?)'
    matches = re.findall(pattern, chemical_formula)
    
    elements = {}
    for match in matches:
        element = match[0]
        quantity = match[1]
        quantity = float(quantity) if quantity else 1.0
        elements[element] = quantity

    return elements

@app.route('/')
def index():
    return render_template('index#2.html')

@app.route('/calculate', methods=['GET','POST'])
def calculate():
    cols = plotly.colors.DEFAULT_PLOTLY_COLORS
    data = request.form
    chemical_formula = data['chemical_formula']
    elements = get_elements(chemical_formula)
    if data['type_energy'] == 'energy':
        energy = float(data['wavelength'])
        energy = energy*1000
    else:
        wavelength = float(data['wavelength'])
        energy = (h*c)/(wavelength*(1e-10))
    cappilary_diameter = data['capillary_diameter']
    # Same calculation logic as in your PyQt app
    try:
        total_mass = sum([elements[element] * xr.atomic_mass(element) for element in elements])
    except ValueError:
        flash('Invalid chemical formula', 'error')
        return redirect('/')
    total_volume = sum([elements[element] * (1e-23) for element in elements])
    packing_fraction = float(data['packing_fraction'])
    density = (total_mass*mu)/total_volume #g/cm³
    packing_density = density*packing_fraction
    distance = float(cappilary_diameter)*(0.1)
    m_u_t = 0
    for element in elements:
        mass_percentage = (elements[element]*xr.atomic_mass(element))/total_mass
        m_u_t += ((mass_percentage)*xr.mu_elam(element, energy)) #(cm²/g)

    m_u_t = m_u_t*packing_density #(cm⁻¹)
    transmission = (math.exp(-(distance*m_u_t)))*100
    mu_R = m_u_t*(distance/2)

    fig = make_subplots(rows=1, cols=2, subplot_titles=(r'$\text{Mass Attenuation Coefficient} \ \frac{\mu}{\rho}$', r'$\mu R \ \text{(Attenuation Coefficient} \times \text{Capillary Radius)}$'))
    i=0
    for element in elements:
        energy_range = np.arange(5000, 30000, 10)
        mu_values = xr.mu_elam(element, energy_range)
        mu_R_values = xr.mu_elam(element, energy_range)*(distance/2)*(density)
        fig.add_trace(go.Scatter(x=energy_range/1000, y=mu_values, line=dict(width=2, color=cols[i]), name=element, showlegend=False), row=1,col=1)
        fig.add_trace(go.Scatter(x=energy_range/1000, y=mu_R_values, line=dict(width=2, color=cols[i]), name=element), row=1, col=2)
        i+=1

    fig.add_scatter(x=[energy/1000], y=[mu_R], name="Sample's µR",row=1, col=2)
    fig.add_hline(y=5, line_dash="dash", line_color ='black', name='µR = 5',row =1, col=2, showlegend=True)
    fig.add_hline(y=1, line_dash="dash", line_color ='blue', name='µR = 1',row=1, col=2, showlegend=True)
            

    fig.update_xaxes(title_text="Energy (keV)", type="log", row=1, col=1)
    fig.update_xaxes(title_text="Energy (keV)", type="log", row=1, col=2)
    fig.update_yaxes(title_text=r"$\mu/\rho \ \text{(cm²/g)}$", type="log", row=1, col=1)
    fig.update_yaxes(title_text=r"$\mu R$", type="log", row=1, col=2)

    density = f'{density:.4f}'
    transmission = f'{transmission:.4f}'
    packing_density = f'{packing_density:.4f}'
    # Return results to be rendered in the HTML
    return jsonify({
        'density': density,
        'packing_density': packing_density, 
        'transmission': transmission,  
        'plot': fig.to_html()
    })

if __name__ == '__main__':
    app.run(debug=True)
