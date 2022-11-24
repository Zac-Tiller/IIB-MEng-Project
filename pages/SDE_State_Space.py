import streamlit as st
from sampling_functions import *

header = st.container()

parameter_selector = st.container()


AMatrixSelector = st.container()

HMatrixSelector = st.container()

theta_selector = st.container()

gamma_simulation_engine = st.container()

NormalGamma_p_selector = st.container()

normal_gamma_simulaton_engine = st.container()

st.markdown(
    """
    <style>
    .main{
    background-color: #F5F5F5;
    }
    <style>
    """,
    unsafe_allow_html=True
)


def RunStateSpaceSimulation(alpha, beta, samples, mean, var, theta):
    gamma_obj = GammaDistr(alpha=alpha, beta=beta)
    gamma_obj.set_process_conditions(t0=0, T=1, sample_size=samples)

    gamma_sim = DistributionSimulator(gamma_obj)

    path, time, jump_time_set = gamma_sim.process_simulation(
        gamma_obj)  # run the simulation. after this comment, the gamma_sim object has the jump_time_sets (sorted process set) as an attribute

    normal_obj = NormalDistr(mean=mean, std=np.sqrt(var), secondary_distr=gamma_obj)
    # normal_gamma_sim = DistributionSimulator(normal_obj)

    state_simulation = StateSpaceSimulation(gamma_sim)

    state_simulation.define_model_h([0, 1])
    state_simulation.define_A_matrix([0, 1, 0, theta])
    state_simulation.set_NG_conditions(normal_obj.mean, (normal_obj.std) ** 2)

    fig1, fig2 = state_simulation.run_state_space_simulation()

    X0_Plot, X1_Plot = st.columns(2)
    X0_Plot.pyplot(fig1)
    X1_Plot.pyplot(fig2)
    return None


with header:
    st.title('Interactive SDE Forward Simulation')

with AMatrixSelector:
    st.markdown('**Choose A Matrix:**')
    st.markdown('-- > Empty for Now. Feature to be built later')


with HMatrixSelector:
    st.markdown('**Choose h Vector:**')
    st.markdown('-- > Empty for Now. Feature to be built later')


with theta_selector:
    st.markdown('**Select Decay Strength in A Matrix**')

    theta = st.slider('Theta Value?', min_value=float(-3), max_value=float(3), step= 0.01)

with parameter_selector:
    st.markdown('**Choose parameters of the Gamma Distribution, Alpha and Beta:**')

    alpha_col, beta_col, sample_col = st.columns(3)

    alpha = alpha_col.slider('Alpha Value?', min_value=float(0.01), max_value=float(5), step=0.01)
    beta = beta_col.slider('Beta Value?', min_value=float(0.01), max_value=float(5), step=0.01)
    samples = sample_col.slider('Number of Samples of The Exp. Distr?', min_value=1, max_value=1000)

with NormalGamma_p_selector:
    st.markdown('**Choose Extra parameters for the Normal-Gamma, Mean and Var:**')

    mean_col, var_col = st.columns(2)

    mean = mean_col.slider('Mean Value?', min_value=float(0.0), max_value=float(5), step=0.01)
    var = var_col.slider('Variance Value?', min_value=float(0.0), max_value=float(5), step=0.01)


with normal_gamma_simulaton_engine:
    if st.button('Click for SDE Forward Simulation'):
        RunStateSpaceSimulation(alpha, beta, samples, mean, var, theta)
