import json
import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from dash import Dash, dcc, html
import plotly.express as px
import plotly.graph_objects as go
from dash.dependencies import Input, Output, State
from plotly.subplots import make_subplots
import io
import base64

TARGET_PTM = {
    'Citrullination': {'name': 'Citrullination',
                       'path': 'data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_citrullination_report-lib.parquet',
                       'symbol': r'UniMod:7', 'site': 'R', 'mass_shift': 0.984016},
    'Hypusine': {'name': 'Hypusine',
                 'path': 'data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_hypusine_report-lib.parquet',
                 'symbol': r'Hypusine', 'site': 'K', 'mass_shift': 114.042927},
    'Deoxyhypusine': {'name': 'Deoxyhypusine',
                      'path': 'data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_deoxyhypusine_report-lib.parquet',
                      'symbol': r'Deoxyhypusine', 'site': 'K', 'mass_shift': 98.031300},
}

# Initialize the Dash app
app = Dash(__name__)

# Define the layout
app.layout = html.Div([
    dcc.Dropdown(
        id='file-dropdown',
        options=[{'label': key, 'value': key} for key in TARGET_PTM.keys()],
        value='Citrullination',
        placeholder='Select File'
    ),

    dcc.Graph(id='matplotlib-plot', style={'height': '800px', 'width': '100%'}),
    dcc.Dropdown(
        id='cancer-type-dropdown',
        options=[],
        value=None,
        multi=True,
        placeholder='Select Cancer Types'
    ),
    dcc.Slider(
        id='bins-slider',
        min=10,
        max=100,
        step=1,
        value=30,
        marks={i: str(i) for i in range(10, 101, 10)},
        tooltip={"placement": "bottom", "always_visible": True}
    ),
    dcc.Graph(id='histogram-plot', style={'height': '800px', 'width': '100%'}),
    dcc.Store(id='processed-data')
])

@app.callback(
    [Output('cancer-type-dropdown', 'options'),
     Output('processed-data', 'data')],
    [Input('file-dropdown', 'value')]
)
def update_dropdowns(selected_ptm):
    if not selected_ptm:
        return [[]], None
    ptm_info = TARGET_PTM[selected_ptm]
    df = pd.read_parquet(ptm_info['path'])
    if df is None:
        return [[]], None
    metadata = pd.read_csv("data/metadata.csv")
    dat = df.merge(metadata, on='ID', how='left')
    dat['has_target_PTM'] = dat['Modified.Sequence'].str.contains(ptm_info['symbol'])
    dat['has_target_site'] = dat['Stripped.Sequence'].str.contains(ptm_info['site'])
    dat['num_PTM'] = dat['Modified.Sequence'].str.count(ptm_info['symbol'])
    total_intensity = dat.groupby(['ID', 'assay', 'has_target_PTM']).agg(total_intensity=('Precursor.Quantity', 'sum')).reset_index()
    total_intensity['total_intensity'] = np.log10(total_intensity['total_intensity'])
    total_intensity = total_intensity.merge(metadata, on='ID')
    cancer_types = [{'label': ct, 'value': ct} for ct in total_intensity['Cancer Type'].dropna().unique()]
    return cancer_types, total_intensity.to_dict('records')

@app.callback(
    Output('histogram-plot', 'figure'),
    [Input('cancer-type-dropdown', 'value'),
     Input('bins-slider', 'value')],
    [State('file-dropdown', 'value'),
     State('processed-data', 'data')]
)
def update_histogram(selected_cancer_types, bins, selected_ptm, processed_data):
    if not selected_cancer_types or not processed_data:
        return go.Figure()
    total_intensity = pd.DataFrame(processed_data)
    ptm_info = TARGET_PTM[selected_ptm]
    total_intensity_labeled = total_intensity[total_intensity['has_target_PTM']].dropna(subset=['Cancer Type'])
    total_intensity_labeled = total_intensity_labeled[total_intensity_labeled['Cancer Type'].isin(selected_cancer_types)]
    total_intensity_labeled = total_intensity_labeled[total_intensity_labeled['assay'].isin(['FT', 'IgB'])]
    unique_assays = total_intensity_labeled['assay'].unique()
    fig = make_subplots(rows=2, cols=1, subplot_titles=unique_assays)
    cmap = plt.get_cmap('viridis')
    unique_cancer_types = total_intensity_labeled['Cancer Type'].unique()
    colors = {ct: matplotlib.colors.rgb2hex(cmap(i / len(unique_cancer_types))) for i, ct in enumerate(unique_cancer_types)}
    for i, assay in enumerate(unique_assays):
        subset_assay = total_intensity_labeled[total_intensity_labeled['assay'] == assay]
        for ct in unique_cancer_types:
            subset = subset_assay[subset_assay['Cancer Type'] == ct]
            fig.add_trace(go.Histogram(
                x=subset['total_intensity'],
                name=ct,
                marker_color=colors[ct],
                opacity=0.5,
                nbinsx=bins,
                showlegend=(i == 0)  # Show legend only for the first subplot
            ), row=(i + 1), col=1)
    fig.update_layout(
        barmode='overlay',
        xaxis_title='Total intensity (log10)',
        yaxis_title='Frequency',
        title=f"Distribution of total intensity of modified peptide with {ptm_info['name']} on site {ptm_info['site']} by Cancer Type",
        legend_title_text='Cancer Type'
    )
    return fig
@app.callback(
    Output('matplotlib-plot', 'figure'),
    [Input('file-dropdown', 'value'),
     Input('processed-data', 'data')]
)
def update_matplotlib_plot(selected_ptm, processed_data):
    if not processed_data:
        return {'data': [], 'layout': {}}

    total_intensity = pd.DataFrame(processed_data)
    ptm_info = TARGET_PTM[selected_ptm]
    total_intensity_labeled = total_intensity.dropna(subset=['group'])
    color_mapping = {
        ('Case', True): 'green',
        ('Case', False): 'orange',
        ('Control', True): 'red',
        ('Control', False): 'yellow'
    }

    unique_assays = total_intensity_labeled['assay'].unique()
    fig = make_subplots(rows=2, cols=1, subplot_titles=unique_assays)

    for i, assay in enumerate(unique_assays):
        subset_assay = total_intensity_labeled[total_intensity_labeled['assay'] == assay]
        for (case_label, mod_label), color in color_mapping.items():
            subset = subset_assay[(subset_assay['group'] == case_label) &
                                  (subset_assay['has_target_PTM'] == mod_label)]
            label = f'{case_label}, {"Modified peptide" if mod_label else "Unmodified peptide"}'
            fig.add_trace(go.Histogram(
                x=subset['total_intensity'],
                name=label,
                marker_color=color,
                opacity=0.5,
                showlegend=(i == 0)  # Show legend only for the first subplot
            ), row=(i + 1), col=1)

    fig.update_layout(
        barmode='overlay',
        xaxis_title='Total intensity (log10)',
        yaxis_title='Frequency',
        title=f'Distribution of total intensity of peptide with {ptm_info["name"]} on site {ptm_info["site"]}',
        legend_title_text='Group'
    )

    return fig

# Run the app
if __name__ == '__main__':
    app.run_server(debug=True)