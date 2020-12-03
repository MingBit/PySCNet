from __future__ import absolute_import
from pathlib import Path
import io
import os
import uuid
import dash_uploader as du
import dash
import flask
import random
from umap import UMAP
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import pandas as pd
from textwrap import dedent
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
import numpy as np
import plotly.express as px
import dash_cytoscape as cyto
import pyscnet.NetEnrich as ne
import pyscnet.Preprocessing as pp
from dash.dependencies import Input, Output

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.LUMEN])
app.title = 'PySCNet Dashboard'
UPLOAD_FOLDER_ROOT = r"/home/mwu/dash-sample-apps/apps/dash-pyscnet/data/"
du.configure_upload(app, UPLOAD_FOLDER_ROOT)

FONT_STYLE = {
    "color": '#fff1e6',
    'font-size': '50',
    'margin-left': '1rem',
    'margin-right': '1rem'
}

SIDEBAR_STYLE = {
    "height": "100%",
    "font-size": "45px",
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "20%",
    "padding": "2rem 1rem",
    "background-color": "#1d3557",
}

# the styles for the main content position it to the right of the sidebar and
# add some padding.
CONTENT_STYLE = {
    "top": 0,
    "bottom": 0,
    "margin-left": "20%",
    "padding": "2rem 1rem",
    # "background-color": "#264653",
    "height": "100%",
    "width": "100%",
    "font-size": '30px'

}

DROPDOWN_STYLE = {
    'height': '30px'
}


def get_upload_component(id):
    return du.Upload(
        id=id,
        max_file_size=1800,  # 1800 Mb
        filetypes=['pk'],
        upload_id=uuid.uuid1(),  # Unique session id
    )


def get_GRN_method():
    assert len(
        list(filter(lambda x: 'links' in x, list(object.NetAttrs.keys())))) != 0, "No valid link table available!"
    return list({'label': i, 'value': i} for i in list(filter(lambda x: 'links' in x, list(object.NetAttrs.keys()))))


def get_gene_rank():
    assert 'centralities' in object.NetAttrs.keys(), "No node centralities available!"
    return list({'label': i, 'value': i} for i in list(object.NetAttrs['centralities'].columns[1:]) + ['None'])


def get_cellinfo_column():
    assert 'CellInfo' in object.CellAttrs.keys(), 'No CellInfo available!'
    return list({'label': i, 'value': i} for i in list(object.CellAttrs['CellInfo'].columns))


def windown_sliding_corr(genes, cell_order=None, cell_filter=None, pairwise=True):
    r_window_size = 100
    df = object.ExpMatrix.T

    if cell_order is not None:
        if cell_filter is not None:
            cells = object.CellAttrs['CellInfo'].sort_values(cell_order, ascending=False)
            df = df.loc[cells[cells[cell_order] == cell_filter].index]
        else:
            df = df.reindex(object.CellAttrs['CellInfo'].sort_values(cell_order, ascending=False).index)

    # Interpolate missing data.
    df_interpolated = df.interpolate()
    # Compute rolling window synchrony
    tmp = pd.DataFrame(df[genes].rolling(window=r_window_size, center=True).mean())

    if pairwise:

        rolling_r = pd.DataFrame(
            df_interpolated[genes[0]].rolling(window=r_window_size, center=True).corr(df_interpolated[genes[1]]))
        rolling_r.columns = ['Correlation']
        return rolling_r, tmp
    else:

        return tmp


def update_object(grn_method, top_links, resolution=0.5):
    new_object = ne.buildnet(object, key_links=grn_method, top=int(top_links))
    new_object = ne.get_centrality(new_object)
    new_object = ne.detect_community(new_object, resolution=resolution)

    return new_object


def update_filter_link(grn_method, top_links, node_size, resolution=0.5):
    object = update_object(grn_method, top_links, resolution)
    filtered_link = object.NetAttrs[grn_method].sort_values('weight', ascending=False).head(top_links)
    gene_module = object.NetAttrs['communities']

    color = ["#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(len(np.unique(gene_module.group)))]

    gene_module['color'] = [color[i] for i in gene_module.group]
    gene_module = pd.concat([object.NetAttrs['centralities'], gene_module[['color']].reset_index(drop=True)],
                            axis=1)
    object.NetAttrs['communities'] = gene_module
    gene_module['size'] = np.repeat(30, gene_module.shape[0]) if node_size == 'None' else gene_module[
                                                                                              node_size] * 100 + 1

    nodes = [{'data': {'id': name, 'label': name, 'betweenness': betweenness, 'closeness': closeness,
                       'degree': degree, 'pageRank': pageRank, 'color': color, 'size': size}} for
             name, betweenness, closeness, degree, pageRank, color, size in
             list(gene_module.itertuples(index=False, name=None))]
    edges = [{'data': {'source': source, 'target': target, 'weight': weight}} for source, target, weight in
             list(filtered_link.itertuples(index=False, name=None))]

    new_elements = nodes + edges

    return new_elements


def update_sub_network(click_node=None):
    click_node = list(object.NetAttrs['graph'].node)[0] if click_node is None else click_node
    neighbours = list(object.NetAttrs['graph'].neighbors(click_node))
    gene_module = object.NetAttrs['communities']

    sub_link = pd.DataFrame({'source': np.repeat(click_node, len(neighbours)),
                             'target': neighbours})
    sub_gene_module = gene_module[gene_module.node.isin([click_node] + neighbours)][['node', 'color']]

    sub_nodes = [{'data': {'id': name, 'label': name, 'color': color}} for name, color in
                 sub_gene_module.itertuples(index=False, name=None)]
    sub_edges = [{'data': {'source': source, 'target': target}} for source, target in
                 list(sub_link.itertuples(index=False, name=None))]

    sub_element = sub_nodes + sub_edges

    color = gene_module.loc[gene_module.node == click_node, 'color'].to_list()
    sub_module = gene_module[gene_module.color == color[0]][['node', 'color']]

    inter_node = set(neighbours) & set(sub_module)
    sub_nodes_2 = [{'data': {'id': name, 'label': name, 'color': color}} for name, color in
                   sub_module.itertuples(index=False, name=None)]
    sub_edges_2 = [{'data': {'source': source, 'target': target}}
                   for source, target in
                   list(sub_link.loc[sub_link.target.isin(inter_node)].itertuples(index=False, name=None))]

    sub_element_module = sub_nodes_2 + sub_edges_2

    return [[click_node] + neighbours, sub_element, sub_element_module]


def update_sub_page(sub_page):
    if 'object' in globals():
        if sub_page == 'page_1':
            new_sub_page = create_page_1()
        elif sub_page == 'page_3':
            new_sub_page = create_page_3()
        else:
            new_sub_page = create_page_2()

    else:
        new_sub_page = html.Div(id='page_1-1',
                                children=html.H1('please upload object!',
                                                 style={'font-size': '80px', 'font-weight': 'bold',
                                                        'color': '#457b9d'}))

    return new_sub_page


def get_gene_curve(df, title, yaxis_title):
    f = px.scatter(df, title=title)

    f.update_layout(xaxis_title="Cell name",
                    yaxis_title=yaxis_title,
                    plot_bgcolor='#264653',
                    font=dict(family="Courier New, monospace",
                              size=50,
                              color="RebeccaPurple"),
                    legend=dict(font=dict(family="Courier", size=50),
                                itemsizing='constant'), legend_title_text='')

    f.update_layout(legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=0.01
    ))
    f.update_xaxes(showgrid=False, showticklabels=False)
    f.update_yaxes(showgrid=False, tickfont=dict(family='Arial', color='black', size=28))

    return f


def create_page_1():
    page_1 = html.Div([
        html.H2('Cell distributed by UMAP:', style={'font-size': '50px', 'font-weight': 'bold', 'color': '#fff1e6'}),
        dbc.Row([
            dbc.Col([

                html.P('set PCA components:', style=FONT_STYLE),
                dcc.Input(id='components', type='number', value=10),
            ]),

            dbc.Col([
                html.P('set K neighbors:', style=FONT_STYLE),
                dcc.Input(id='neighbors', type='number', value=300),
            ]),

            dbc.Col([

                html.P('choose color 1', style=FONT_STYLE),
                dcc.Dropdown(options=get_cellinfo_column(), id='color_code_1',
                             value=object.CellAttrs['CellInfo'].columns[0],
                             style=DROPDOWN_STYLE)
            ]),

            dbc.Col([
                html.P('choose color 2', style=FONT_STYLE),
                dcc.Dropdown(options=get_cellinfo_column(), id='color_code_2',
                             value=object.CellAttrs['CellInfo'].columns[1],
                             style=DROPDOWN_STYLE)
            ]),

            dbc.Col([
                html.P('Choose 2D or 3D:', style=FONT_STYLE),
                dcc.RadioItems(id='xD',
                               options=list({'label': i, 'value': i} for i in ['2D', '3D']),
                               value='3D', labelStyle={'display': 'inline-block', 'margin-right': '1rem'},
                               style=FONT_STYLE)
            ])
        ], style={'height': '200px', 'margin-left': '1rem', 'margin-right': '1rem'}),

        dbc.Row([
            dbc.Col([
                html.H2('Cell distribution encoded by color 1',
                        style={'font-size': '50px', 'font-weight': 'bold', 'color': '#fff1e6'}),
                html.Br(),
                dcc.Graph(id='cell_distribution_1', style={'height': '1500px'})
            ]),
            dbc.Col([html.H2('Cell distribution encoded by color 2',
                             style={'font-size': '50px', 'font-weight': 'bold', 'color': '#fff1e6'}),
                     html.Br(),
                     dcc.Graph(id='cell_distribution_2', style={'height': '1500px'})
                     ])
        ], style={'margin-left': '1rem', 'margin-right': '1rem'}),

        dbc.Row([
            html.A('download cell annotation table', style={'margin-left': '3rem', 'color': '#fff1e6',
                                                            'font-size': '60'},
                   id='cell_download', href='/home/CellInfo-download.xlsx')
        ])

    ],
        id='wrapper', style={"background-color": "#264653"})
    return page_1


def create_page_2():
    page_2 = html.Div([
        html.H2('Choose two genes:', style={'font-size': '50px', 'font-weight': 'bold', 'color': '#fff1e6'}),
        dbc.Row([
            dbc.Col([
                html.P('Gene A', style=FONT_STYLE),
                dcc.Dropdown(options=list({'label': i, 'value': i} for i in list(object.GeneAttrs['GeneInfo'].index)),
                             id='gene_a', value=list(object.GeneAttrs['GeneInfo'].index)[0],
                             style=DROPDOWN_STYLE)
            ]),

            dbc.Col([
                html.P('Gene B', style=FONT_STYLE),
                dcc.Dropdown(options=list({'label': i, 'value': i} for i in list(object.GeneAttrs['GeneInfo'].index)),
                             id='gene_b', value=list(object.GeneAttrs['GeneInfo'].index)[1],
                             style=DROPDOWN_STYLE)
            ]),

            dbc.Col([
                html.P('Cells ordered by ', style=FONT_STYLE),
                dcc.Dropdown(options=get_cellinfo_column(), id='cell_order')
            ]),

            dbc.Col([
                html.P('Cells selected by', style=FONT_STYLE),
                dcc.Dropdown(id='cell_filter')
            ])
        ], style={'height': '200px', 'margin-left': '1rem', 'margin-right': '1rem'}),
        html.Br(),
        dbc.Row([
            dbc.Col([dcc.Graph(id='gene_correlation_1', style={'height': '1500px'})]),
            dbc.Col([dcc.Graph(id='gene_correlation_2', style={'height': '1500px'})]),
        ], style={'margin-left': '1rem', 'margin-right': '1rem'}),
        html.Br(),
        dbc.Row([
            html.A('download gene annotation table', href='/home/GeneInfo-download.xlsx',
                   style={'margin-left': '3rem', 'color': '#fff1e6',
                          'font-size': '60'}, id='gene_download')
        ])
    ], id='wrapper', style={"background-color": "#264653"})

    return page_2


def create_page_3():
    elements = update_filter_link(grn_method=get_GRN_method()[0]['value'], top_links=50,
                                  node_size='None')
    neighbours, sub_element_1, sub_element_2 = update_sub_network(click_node=None)
    rolling_neighbour = windown_sliding_corr(neighbours, pairwise=False)
    gene_dynamics = get_gene_curve(rolling_neighbour, title='Gene expression level', yaxis_title='Expression level')

    def_text = 'please click on the gene node!'

    def_stylesheet = [{
        'selector': 'node',
        'style': {'label': 'data(id)',
                  'color': '#fff1e6',
                  'background-color': 'data(color)',
                  'width': 'data(size)',
                  'height': 'data(size)'
                  }
    }]
    page_3 = html.Div([
        html.H1("Create you own gene network",
                style={'font-size': '100px', 'font-weight': 'bold'}),
        html.Br(),
        dbc.Row([
            dbc.Col([
                cyto.Cytoscape(
                    id='gene_network',
                    layout={'name': 'cose'},
                    style={'width': '95%', 'height': '1500px',
                           'background-color': '#264653'},
                    stylesheet=def_stylesheet,
                    elements=elements
                )], style={'padding': '10px'}),

            dbc.Col([
                html.P('Choose GRN method', style=FONT_STYLE),
                dcc.Dropdown(id='grn_method', options=get_GRN_method(),
                             value=get_GRN_method()[0]['value'], style=FONT_STYLE),

                html.Br(),
                html.P('Choose top links', style=FONT_STYLE),
                dcc.RadioItems(id='top_links',
                               options=list({'label': str(i), 'value': str(i)} for i in [50, 100, 200, 500, 1000]),
                               value="50", labelStyle={'display': 'inline-block', 'margin-right': '1rem'},
                               style=FONT_STYLE),

                html.Br(),
                html.P('Choose resolution for module detection', style=FONT_STYLE),
                dcc.Slider(id='resolution', min=0, max=1, step=0.1, value=0.5,
                           marks={0: {'label': '0'},
                                  0.5: {'label': '0.5'},
                                  1: {'label': 1}}),

                html.Br(),
                html.P('Choose network layout', style=FONT_STYLE),
                dcc.Dropdown(
                    id='net_layout',
                    value='cose',
                    clearable=False,
                    options=[
                        {'label': name.capitalize(), 'value': name}
                        for name in ['cose', 'random', 'circle', 'grid', 'concentric']
                    ], style=FONT_STYLE),

                html.Br(),
                html.P('Node size encoded by', style=FONT_STYLE),
                dcc.Dropdown(id='node_size_encode', options=get_gene_rank(),
                             value='None', style=FONT_STYLE),
                html.Br(),
                html.A('download network table', href='/home/NetInfo-download.xlsx',
                       style={'color': '#fff1e6', 'margin-left': '1rem',
                              'font-size': '60'}, id='net_download')

            ], width=3.8, style={'background-color': '#264653', 'margin-right': '1rem'})

        ], style={"margin-left": "1rem"}),
        html.Br(),
        html.Hr(),
        dbc.Row([
            dbc.Col([
                html.H2(id='node_neighbors', children=def_text,
                        style={'font-size': '50px', 'font-weight': 'bold'}),
                cyto.Cytoscape(
                    id='selected_node_neighbors',
                    layout={'name': 'cose'},
                    style={'width': '92%', 'height': '1000px', 'background-color': '#264653'},
                    stylesheet=def_stylesheet,
                    elements=sub_element_1
                )]),
            dbc.Col([
                html.H2(id='node_module', children=def_text,
                        style={'font-size': '50px', 'font-weight': 'bold'}),
                cyto.Cytoscape(
                    id='selected_node_module',
                    layout={'name': 'grid'},
                    style={'width': '92%', 'height': '1000px', 'background-color': '#264653'},
                    stylesheet=def_stylesheet,
                    elements=sub_element_2
                )
            ])
        ], style={"margin-left": "1rem"}),
        html.Br(),
        html.Hr(),
        html.H2('Gene kinetics',
                style={'font-size': '50px', 'font-weight': 'bold'}),
        dcc.Graph(id='gene_dynamics_plot', figure=gene_dynamics,
                  style={'height': '1500px', 'margin-left': '3rem'})

    ], style={"margin-left": "1rem"}, id='wrapper')

    return page_3


sidebar = html.Div(
    [
        html.H3("PySCNet Dashboard", className="display-4",
                style={'color': 'white', "font-size": "80px", 'font-weight': 'bold'}),
        html.Br(),
        html.P("A python dashboard for pyscnet visualization.",
               style={'color': 'white', "font-size": "30px", 'font-style': 'italic'}),
        html.Hr(),
        dbc.Nav(
            [
                dbc.NavLink("Introduction", href="/page-0", id="page-0-link",
                            style={"margin-bottom": "5rem"}),
                dbc.NavLink("Cell Attributes", href="/page-1", id="page-1-link",
                            style={"margin-bottom": "5rem"}),
                dbc.NavLink("Gene Attributes", href="/page-2", id="page-2-link",
                            style={"margin-bottom": "5rem"}),
                dbc.NavLink("Network Attributes", href="/page-3", id="page-3-link",
                            style={"margin-bottom": "5rem"}),
                dbc.NavLink("Contact", href="/page-4", id="page-4-link",
                            style={"margin-bottom": "5rem"}),
                html.Hr(),
                html.Div(
                    [
                        get_upload_component(id='dash-uploader'),
                        html.H5(id='callback-output', style={'color': 'white'}),
                    ],
                    style={  # wrapper div style
                        'textAlign': 'center',
                        'width': '550px',
                        'padding': '5px',
                        'margin-left': '1rem',
                        'display': 'inline-block',
                        'color': 'white'
                    })
            ],
            vertical=True,
            pills=True
        ),
    ],

    style=SIDEBAR_STYLE,
)

content = html.Div(id="page-content", style=CONTENT_STYLE)
app.layout = html.Div([dcc.Location(id="url"), sidebar, content], style={"display": "flex"})


@app.callback(
    [Output(f"page-{i}-link", "active") for i in range(0, 5)],
    [Input("url", "pathname")])
def toggle_active_links(pathname):
    if pathname == "/":
        return True, False, False, False, False
    return [pathname == f"/page-{i}" for i in range(0, 5)]


@app.callback(Output("page-content", "children"), [Input("url", "pathname")])
def render_page_content(pathname):
    if pathname in ["/", "/page-0"]:
        return dcc.Markdown(dedent(open('assets/Introduction.md', 'r').read()), style={'font-size': '35px'})
    elif pathname == "/page-1":
        return update_sub_page('page_1')
    elif pathname == "/page-2":
        return update_sub_page('page_2')
    elif pathname == "/page-3":
        return update_sub_page('page_3')
    elif pathname == "/page-4":
        return dcc.Markdown(dedent(open('assets/Contact.md', 'r').read()))
    return dbc.Jumbotron(
        [
            html.H1("404: Not found", className="text-danger"),
            html.Hr(),
            html.P(f"The pathname {pathname} was not recognised..."),
        ]
    )


@app.callback(Output('callback-output', 'children'),
              [Input('dash-uploader', 'isCompleted'),
               Input('dash-uploader', 'fileNames'),
               Input('dash-uploader', 'upload_id')])
def get_a_list(is_completed, filenames, upload_id):
    if is_completed and filenames is not None:
        global object
        object = pp.load_Gnetdata_object(UPLOAD_FOLDER_ROOT + '/' + str(upload_id) + '/' + filenames[0])
        return [filenames]


@app.callback(Output('cell_filter', 'options'),
              Input('cell_order', 'value'))
def update_cell_filter(cell_order):
    return list({'label': i, 'value': i} for i in list(np.unique(object.CellAttrs['CellInfo'][cell_order])))


@app.callback(Output('cell_filter', 'value'),
              Input('cell_filter', 'options'))
def initialize_cell_filter(cell_filter_options):
    # return cell_filter_options[0]['value']
    return None


@app.callback([Output('cell_distribution_1', 'figure'),
               Output('cell_distribution_2', 'figure')],
              [Input('xD', 'value'),
               Input('color_code_1', 'value'),
               Input('color_code_2', 'value'),
               Input('components', 'value'),
               Input('neighbors', 'value')])
def update_cell_distribution(xD, color_1, color_2, components, neighbours):
    pca = PCA(n_components=components).fit_transform(object.ExpMatrix.T)
    if xD == '2D':
        proj_2d = UMAP(n_neighbors=neighbours,
                       n_components=2,
                       min_dist=0.8,
                       metric='correlation').fit_transform(pca)
        # proj_2d = TSNE(n_components=2).fit_transform(pca)
        fig_2d_1 = px.scatter(
            proj_2d, x=0, y=1,
            color=object.CellAttrs['CellInfo'][color_1], labels={'color': color_1})

        fig_2d_1.update_layout(legend=dict(font=dict(family="Courier", size=50), itemsizing='constant'),
                               legend_title_text='')
        fig_2d_2 = px.scatter(
            proj_2d, x=0, y=1,
            color=object.CellAttrs['CellInfo'][color_2], labels={'color': color_2})

        fig_2d_2.update_layout(legend=dict(font=dict(family="Courier", size=50), itemsizing='constant'),
                               legend_title_text='')
        return fig_2d_1, fig_2d_2
    else:
        proj_3d = UMAP(n_neighbors=neighbours,
                       n_components=3,
                       min_dist=0.5,
                       metric='correlation').fit_transform(pca)

        fig_3d_1 = px.scatter_3d(
            proj_3d, x=0, y=1, z=2,
            color=object.CellAttrs['CellInfo'][color_1], labels={'color': color_1})

        fig_3d_1.update_layout(legend=dict(font=dict(family="Courier", size=50), itemsizing='constant'),
                               legend_title_text='')

        fig_3d_2 = px.scatter_3d(
            proj_3d, x=0, y=1, z=2,
            color=object.CellAttrs['CellInfo'][color_2], labels={'color': color_2})
        fig_3d_2.update_layout(legend=dict(font=dict(family="Courier", size=50), itemsizing='constant'),
                               legend_title_text='')
        return fig_3d_1, fig_3d_2


@app.callback([Output('gene_correlation_1', 'figure'),
               Output('gene_correlation_2', 'figure')],
              [Input('gene_a', 'value'),
               Input('gene_b', 'value'),
               Input('cell_order', 'value'),
               Input('cell_filter', 'value')])
def gene_cor_curve(gene_a, gene_b, cell_order, cell_filter):
    rolling_1, rolling_2 = windown_sliding_corr([gene_a, gene_b], cell_order, cell_filter, pairwise=True)
    f1 = get_gene_curve(rolling_1, title='window_sliding correlation', yaxis_title='Correlation')
    f2 = get_gene_curve(rolling_2, title='Gene expression level', yaxis_title='Expression level')
    return f2, f1


@app.callback([Output('gene_network', 'elements'),
               Output('gene_network', 'layout'),
               Output('gene_network', 'stylesheet')],
              [Input('grn_method', 'value'),
               Input('top_links', 'value'),
               Input('net_layout', 'value'),
               Input('resolution', 'value'),
               Input('node_size_encode', 'value')])
def update_gene_network(grn_method, top_links, net_layout, resolution, node_size_encode):
    new_elements = update_filter_link(grn_method=grn_method,
                                      top_links=int(top_links),
                                      node_size=node_size_encode,
                                      resolution=resolution)

    new_layout = {'name': net_layout, 'animate': True}

    new_stylesheet = [
        {
            'selector': 'node',
            'style': {
                'label': 'data(id)',
                'background-color': 'data(color)',
                'width': 'data(size)',
                'height': 'data(size)',
                'color': '#fff1e6'}
        }]
    return [new_elements, new_layout, new_stylesheet]


@app.callback([Output('selected_node_neighbors', 'elements'),
               Output('selected_node_module', 'elements'),
               Output('selected_node_neighbors', 'stylesheet'),
               Output('selected_node_module', 'stylesheet'),
               Output('node_neighbors', 'children'),
               Output('node_module', 'children'),
               Output('gene_dynamics_plot', 'figure')],
              [Input('gene_network', 'tapNodeData')])
def update_sub_net(data):
    if data:
        neighbours, new_sub_elements_1, new_sub_elements_2 = update_sub_network(data['id'])
        new_stylesheet_1 = [{
            'selector': 'node',
            'style': {
                'label': 'data(id)',
                'color': '#fff1e6',
                'background-color': 'data(color)'
            }
        }]

        neighbour_text = 'Genes connected to ' + data['id']
        module_text = 'Genes assigned to the same module as ' + data['id']

        rolling_neighbour = windown_sliding_corr(neighbours, pairwise=False)
        gene_dynamics = get_gene_curve(rolling_neighbour, title='Gene expression level', yaxis_title='Expression level')
    return [new_sub_elements_1, new_sub_elements_2, new_stylesheet_1,
            new_stylesheet_1, neighbour_text, module_text, gene_dynamics]


@app.server.route('/home/<path:path>')
def download_excel(path):
    download_path = os.getenv('HOME') + '/Desktop/'
    if path in ['GeneInfo-download.xlsx', 'CellInfo-download.xlsx']:

        if path == 'GeneInfo-download.xlsx':
            df = object.GeneAttrs['GeneInfo']
        else:
            df = object.CellAttrs['CellInfo']

        absolute_filename = os.path.join(download_path, path)
        writer = pd.ExcelWriter(absolute_filename)
        df.to_excel(writer, 'Sheet1')
        writer.save()

    elif path == 'NetInfo-download.xlsx':
        absolute_filename = os.path.join(download_path, path)
        writer = pd.ExcelWriter(absolute_filename)
        for name in list(filter(lambda x: 'links' in x, list(object.NetAttrs.keys()))):

            df = object.NetAttrs[name]
            df.to_excel(writer, name)
            writer.save()
    else:
        print('404 error!')

    return flask.send_from_directory(download_path, path)


if __name__ == '__main__':
    app.run_server(host='127.0.0.1', port='8080', debug=True,
                   dev_tools_ui=False, dev_tools_props_check=False)
