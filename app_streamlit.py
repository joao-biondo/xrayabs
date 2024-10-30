import streamlit as st
import xraydb as xr
import numpy as np
import plotly.graph_objects as go
import re
import scipy.constants
import math, plotly
from plotly.subplots import make_subplots
import base64

st.set_page_config(page_title="X-ray Attenuation Calculator", page_icon="Paineira#2.png", layout="wide")
st.title("X-ray Attenuation Calculator")
st.markdown(r'''This page was desinged for X-ray attenuation calculation for samples and experiments performed at the Paineira beamline. 
                The user is only required to insert the sample's chemical formula, the X-ray energy or wavelength and select the capillary diameter and packing fraction.
                Each item has a description that should help filling the form correctly. Also, the calculations are explained at the end of the page.
                **OBS:** The options given for the capillary diameter are the ones available at the Paineira beamline.
            ''')


def get_img_as_base64(file):
    with open(file, "rb") as f:
        data = f.read()
    return base64.b64encode(data).decode()

img = get_img_as_base64("Paineira.png")

# Constantes
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

def calculate(chemical_formula, energy_or_wavelength, type_energy, capillary_diameter, packing_fraction):
    elements = get_elements(chemical_formula)
    if type_energy == 'Energy (keV)':
        energy = float(energy_or_wavelength) * 1000
    else:
        wavelength = float(energy_or_wavelength)
        energy = (h*c)/(wavelength*(1e-10))
    distance = float(capillary_diameter.split(sep=' ')[0])*0.1

    
    total_mass = sum([elements[element] * xr.atomic_mass(element) for element in elements])
    

    total_volume = sum([elements[element] * (1e-23) for element in elements])
    density = (total_mass * mu) / total_volume
    packing_density = density * packing_fraction

    m_u_t = sum(
        ((elements[element] * xr.atomic_mass(element)) / total_mass) * xr.mu_elam(element, energy)
        for element in elements
    ) * packing_density
    transmission = math.exp(-distance * m_u_t) * 100
    mu_R = m_u_t * (distance / 2)

    return density, packing_density, transmission, energy, mu_R, distance, total_mass

def test_chemical_element(chemical_formula):
    elements = get_elements(chemical_formula)

page_bg_img = f"""
<style>
[data-testid="stMain"] {{
    background-image: url("data:image/png;base64,{img}");
    background-size: cover;

}}

[data-testid="stHeader"] {{
    background-color: rgba(256,0,0,0);

}}

</style>
"""

st.markdown(page_bg_img, unsafe_allow_html=True)


# Entradas do Usuário com Validação
chemical_formula = st.text_input("Enter the sample's chemical formula. Be aware that capitalization is required (ex: YBa2Cu3O6.5)")
if chemical_formula and not re.match(r"^([A-Z][a-z]?\d*\.?\d*)+$", chemical_formula):
    st.error("Invalid chemical Formula. Capitalization is required (ex: YBa2Cu3O6.5).")

energy_or_wavelength = st.text_input("Enter the X-ray energy in keV or the wavelength in Å")
if energy_or_wavelength and not re.match(r"^-?\d+(\.\d+)?$", energy_or_wavelength):
    st.error("Enter the X-ray energy in keV or the wavelength in Å.")

type_energy = st.selectbox("Select the type of entry", options=["Wavelength (Å)", "Energy (keV)"])
capillary_diameter = st.selectbox("Capillary Diameter (mm)", [
                "1.00 mm - Kapton",
                "0.30 mm - Kapton",
                "0.50 mm - Kapton",
                "0.70 mm - Kapton",
                "1.12 mm - Kapton",
                "1.37 mm - Kapton",
                "1.57 mm - Kapton",
                "0.80 mm (ID) - Quartzo",
                "1.00 mm (ID) - Quartzo",
                "1.50 mm (ID) - Quartzo",
                "2.00 mm (ID) - Quartzo"
])
packing_fraction = st.text_input("Enter the Packing Fraction.This value represents the decrease in the sample's density when filling the capillary. It should be a value between 0 and 1, and it is often 0.6.")
if packing_fraction and not re.match(r"^0(\.\d+)?|1$", packing_fraction):
    st.error("It must be a value between 0 and 1.")

# Executar Cálculo ao Clicar no Botão
if st.button("Calculate"):
    if chemical_formula and energy_or_wavelength and packing_fraction:
        elements = get_elements(chemical_formula)
        if True:#try:
            test_chemical_element(chemical_formula)
            density, packing_density, transmission, energy, mu_R, distance, total_mass = calculate(
            chemical_formula, energy_or_wavelength, type_energy, capillary_diameter, float(packing_fraction)
        )
        
            st.write(f"Density: {density:.4f} g/cm³")
            st.write(f"Packed Density: {packing_density:.4f} g/cm³")
            st.write(f"µR: {mu_R:.4f}")
            st.write(f"Transmission: {transmission:.2f} %")
            st.write(f"Energy: {energy*(1e-3):.4f} keV")

            # Gráficos

            energy_range = np.arange(5000, 30000, 10)
            mu_list = []
            cols = plotly.colors.DEFAULT_PLOTLY_COLORS
            fig = make_subplots(rows=1, cols=2, subplot_titles=('Mass Attenuation Coefficient - µ/ρ', 'µR (Attenuation Coefficient x Capillary Radius)'))
            i=0
            for element in elements:
                mu_values = xr.mu_elam(element, energy_range)
                mu_list += mu_values*((elements[element] * xr.atomic_mass(element)) / total_mass)
                mu_R_values = xr.mu_elam(element, energy_range)*(distance/2)*(density)
                fig.add_trace(go.Scatter(x=energy_range/1000, y=mu_values, line=dict(width=2, color=cols[i]), name=element, showlegend=False), row=1,col=1)
                fig.add_trace(go.Scatter(x=energy_range/1000, y=mu_R_values, line=dict(width=2, color=cols[i]), name=element), row=1, col=2)
                i+=1

            fig.add_trace(go.Scatter(x=energy_range/1000, y=mu_list, line=dict(width=2, color=cols[i+1]), name="µ Total - Sample", showlegend=False), row=1,col=1)
            fig.add_trace(go.Scatter(x=energy_range/1000, y=(mu_list*(distance/2)*density), line=dict(width=2, color=cols[i+1]), name="µR - Sample"), row=1, col=2)
            
            #fig.add_scatter(x=[energy/1000], y=[mu_R], name="Sample's µR",row=1, col=2)
            fig.add_vline(x=energy/1000, line_dash="dash", line_color ='red', name='µR = 5',row =1, col=2, showlegend=True)
            fig.add_hline(y=5, line_dash="dash", line_color ='black', name='µR = 5',row =1, col=2, showlegend=True)
            fig.add_hline(y=1, line_dash="dash", line_color ='blue', name='µR = 1',row=1, col=2, showlegend=True)
                    
            fig.update_annotations(font=dict(size=25, color='black'))
            fig.update_xaxes(title_font_color='black', title_text="Energy (keV)", type="log", gridcolor='Black', tickfont=dict(color='black'), tickcolor='black', row=1, col=1)
            fig.update_xaxes(title_font_color='black', title_text="Energy (keV)", type="log", gridcolor='Black', tickfont=dict(color='black'), tickcolor='black', row=1, col=2)
            fig.update_yaxes(title_font_color='black', title_text=r"µ/ρ (cm²/g)", type="log", gridcolor='Black', tickfont=dict(color='black'), tickcolor='black', row=1, col=1)
            fig.update_yaxes(title_font_color='black', title_text=r"µR", type="log", gridcolor='Black', tickfont=dict(color='black'), tickcolor='black', row=1, col=2)
            fig.update_layout(legend=dict(title_font_family="Serif", font=dict(size=23)))
            st.plotly_chart(fig)

            # Explicação dos Gráficos
            st.markdown(r"""
            # Graphs and Calculations
            ## Density
            The density is calculated considering that each atom of the unit cell occupies 1 Å³. This is a just an estimate value, but surves our purposes well.
            The Packing Fraction corresponds to the decrease in the sample's density when filling the capillary, so the Packed Density is the product of the Density and the Packing Fraction.
            ## Mass Attenuation Coefficient
            The Mass Attenuation Coefficient of an element gives the probability of interaction of the X-ray photon for a given energy inside a medium composed by that element, per unit of distance traveled.
            The total Mass Attenuation Coefficient is estimated by summing mass attenuation coeficieint of each element multiplied by its mass percentage:
            $$(\frac{\mu}{\rho})_{T} = \sum_{i} (\frac{\mu}{\rho})_{i} \times w_{i}$$, 
                where $w_{i}$ is the mass percentage of the i-element.
            ## Transmission
            The transmission is calculated using the exponential attenuation formula:
            $$\text{Transmission (\%)} = 100 \times e^{-(\frac{\mu}{\rho}) \ \rho \ 2r}$$,
            where μ is the total mass attenuation coefficient, ρ is the sample's density, and r is the capillary radius.
            ## Graphs        
            In the previous section we discussed how to calculate transmission of the X-ray beam given the sample's composition.
            The graphs above show the Mass Attenuation Coeficient (left) and the $\mu R$ (right) value of each element as a function of energy.
            The $\mu R$ value is the product of the mass attenuation coefficient and the capillary radius. This value corresponds to (half of) the argument of the exponential attenuation formula.
            The black dotted line in the right graph represents a $\mu R$ value of 5, whihch gives a transmission of approximately 0.005%. Samples with this kind of attenuation have no sensible X-ray Diffraction signal.
            The blue dotted line represents a $\mu R$ value of 1, which gives a transmission of approximately 13.5%.
            The optimal $\mu R$ value for X-ray Diffraction experiments lies between those two dotted lines.
            """)
        else:#except ValueError:
            st.error("Invalid chemical element")
    else:
        st.warning("Please, fill all the field correctly.")
    
        
