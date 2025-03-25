import streamlit as st
import pandas as pd
import networkx as nx
from pyvis.network import Network
import io
import json
import os

# Function to load predefined datasets
def load_predefined_dataset(dataset_name):
    if dataset_name == "TRRUST Human":
        file_path = "trrust_rawdata.human.tsv"
        if os.path.exists(file_path):
            try:
                # TRRUST format: TF Gene Interaction PubMed_ID
                df = pd.read_csv(file_path, sep="\t", header=None)
                df.columns = ["Source", "Target", "Interaction", "Citation"]
                # Convert interaction types to numeric (1 for activation, 2 for inhibition)
                df["Type"] = df["Interaction"].map({"Activation": 1, "Repression": 2})
                return df
            except Exception as e:
                st.error(f"Error loading TRRUST dataset: {e}")
                return None
        else:
            st.error(f"File not found: {file_path}")
            return None
            
    elif dataset_name == "KEGG Human":
        file_path = "new_kegg.human.reg.direction.txt"
        if os.path.exists(file_path):
            try:
                # Assuming format similar to TRRUST
                df = pd.read_csv(file_path, sep="\t")
                # Make sure it has the required columns
                if "Source" not in df.columns or "Target" not in df.columns:
                    # Rename columns if needed
                    df.columns = ["Source", "Target", "Type", "Citation"] if len(df.columns) >= 4 else ["Source", "Target", "Type"]
                # Ensure Type is numeric
                if "Type" in df.columns:
                    df["Type"] = pd.to_numeric(df["Type"], errors="coerce").fillna(1).astype(int)
                return df
            except Exception as e:
                st.error(f"Error loading KEGG dataset: {e}")
                return None
        else:
            st.error(f"File not found: {file_path}")
            return None
    
    return None

# Function to find interaction between specific genes
def find_gene_interaction(source, target, dataset_name):
    df = load_predefined_dataset(dataset_name)
    if df is None:
        return None
    
    # Look for exact matches
    interactions = df[(df["Source"] == source) & (df["Target"] == target)]
    return interactions if not interactions.empty else None

# Function to draw the network
def draw_network(df, include_add_tools=False):
    G = nx.DiGraph()
    
    # Add nodes from dataframe
    all_nodes = set(df['Source']).union(set(df['Target']))
    for node in all_nodes:
        G.add_node(node, label=node)
    
    # Add edges from dataframe
    for index, row in df.iterrows():
        source = row['Source']
        target = row['Target']
        edge_type = int(row['Type']) if not isinstance(row['Type'], int) else row['Type']
        
        citation = row.get('Citation', "") if 'Citation' in df.columns else ""
            
        if edge_type == 1:
            G.add_edge(source, target, type='activation', citation=citation)
        elif edge_type == 2:
            G.add_edge(source, target, type='inhibition', citation=citation)
    
    # Initialize session state for node/edge management if needed
    if include_add_tools and 'network_nodes' not in st.session_state:
        st.session_state.network_nodes = list(all_nodes)
    
    # Add any additional nodes/edges from session state
    if 'added_nodes' in st.session_state:
        for node in st.session_state.added_nodes:
            if node not in G.nodes():
                G.add_node(node, label=node)
    
    if 'added_edges' in st.session_state:
        for edge in st.session_state.added_edges:
            source, target, edge_type = edge
            if edge_type == 1:
                G.add_edge(source, target, type='activation')
            else:
                G.add_edge(source, target, type='inhibition')
    
    # Create network
    net = Network(notebook=True, height="750px", width="100%", directed=True)
    net.from_nx(G)
    
    # Configure nodes and edges
    for node in net.nodes:
        node['shape'] = 'circle'
        node['size'] = 45
        node['font'] = {'size': 12, 'face': 'arial', 'bold': True, 'color': '#000000'}
        node['color'] = {'border': '#2B7CE9', 'background': '#D2E5FF'}
        node['label'] = node['id']
        node['scaling'] = {'label': {'enabled': False}}
        node['widthConstraint'] = {'minimum': 40, 'maximum': 40}

    for edge in net.edges:
        if edge['type'] == 'activation':
            edge['arrows'] = {'to': {'enabled': True, 'type': 'arrow'}}
            edge['color'] = 'red'
        elif edge['type'] == 'inhibition':
            edge['arrows'] = {'to': {'enabled': True, 'type': 'bar'}}
            edge['color'] = 'blue'
            edge['arrowStrikethrough'] = True
        
        if 'citation' in edge and edge['citation']:
            edge['title'] = edge['citation']
    
    # Configure options
    options_str = """
    {
        "physics": {"enabled": false},
        "nodes": {
            "shape": "circle",
            "size": 45,
            "font": {"size": 12, "face": "arial", "bold": true, "color": "#000000"},
            "scaling": {"label": {"enabled": false}},
            "widthConstraint": {"minimum": 40, "maximum": 40},
            "borderWidth": 2,
            "shadow": true
        },
        "edges": {
            "width": 2,
            "shadow": true,
            "smooth": false
        },
        "interaction": {
            "dragNodes": true,
            "dragView": true,
            "zoomView": true
        }
    }
    """
    
    try:
        options_dict = json.loads(options_str)
        net.set_options(json.dumps(options_dict))
    except Exception as e:
        st.warning(f"Could not set network options: {e}")
    
    return net, G

# Initialize session state
if 'added_nodes' not in st.session_state:
    st.session_state.added_nodes = []
if 'added_edges' not in st.session_state:
    st.session_state.added_edges = []
if 'df' not in st.session_state:
    st.session_state.df = None

# UI
st.title("Gene Regulatory Network Visualization")

tab1, tab2, tab3 = st.tabs(["File Upload", "Direct Input", "Predefined Datasets"])

with tab1:
    uploaded_file = st.file_uploader("Upload TSV file", type=["tsv", "txt"])
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file, sep="\t")
        st.session_state.df = df
        st.dataframe(df.head())

with tab2:
    example = "Source\tTarget\tType\tCitation\ngeneA\tgeneB\t1\tCitation 1\ngeneB\tgeneC\t2\tCitation 2"
    data = st.text_area("Paste TSV data:", value=example, height=200)
    if st.button("Parse Input"):
        try:
            df = pd.read_csv(io.StringIO(data), sep="\t")
            st.session_state.df = df
            st.dataframe(df.head())
        except Exception as e:
            st.error(f"Error parsing data: {e}")

with tab3:
    st.subheader("Predefined Datasets")
    dataset_option = st.selectbox(
        "Select dataset:",
        ["TRRUST Human", "KEGG Human"],
        key="predefined_dataset"
    )
    
    if st.button("Load Full Dataset"):
        df = load_predefined_dataset(dataset_option)
        if df is not None:
            st.session_state.df = df
            st.success(f"Loaded {dataset_option} dataset with {len(df)} interactions")
            st.dataframe(df.head())
    
    # Gene interaction search section
    st.subheader("Search Gene Interaction")
    col1, col2, col3 = st.columns([2, 2, 1])
    
    with col1:
        lookup_source = st.text_input("Source Gene:", key="lookup_source")
    with col2:
        lookup_target = st.text_input("Target Gene:", key="lookup_target")
    with col3:
        lookup_dataset = st.selectbox(
            "Dataset:",
            ["TRRUST Human", "KEGG Human"],
            key="lookup_dataset"
        )
    
    if st.button("Find Interaction", key="find_interaction"):
        if lookup_source and lookup_target:
            interaction = find_gene_interaction(lookup_source, lookup_target, lookup_dataset)
            if interaction is not None:
                st.success(f"Found interaction between {lookup_source} and {lookup_target}")
                st.dataframe(interaction)
                
                # Add to current network
                if st.session_state.df is None:
                    st.session_state.df = interaction
                    st.success("Added to the network.")
                else:
                    # Check if interaction already exists
                    exists = any((st.session_state.df["Source"] == lookup_source) & 
                                 (st.session_state.df["Target"] == lookup_target))
                    
                    if not exists:
                        st.session_state.df = pd.concat([st.session_state.df, interaction], ignore_index=True)
                        st.success("Added to the network.")
                    else:
                        st.info("This interaction already exists in the network.")
                
                st.info("Press 'Generate Network' to update the visualization.")
            else:
                st.warning(f"No interaction found between {lookup_source} and {lookup_target} in {lookup_dataset}.")
        else:
            st.warning("Please enter both source and target genes.").

# Visualization section
if st.session_state.df is not None:
    st.subheader("Network Visualization")
    
    # Main network display
    if st.button("Generate Network", key="generate"):
        net, G = draw_network(st.session_state.df, include_add_tools=True)
        
        # Save the graph structure
        st.session_state.network_nodes = list(G.nodes())
        
        # Display network
        net.save_graph("network.html")
        with open("network.html", 'r', encoding='utf-8') as f:
            html_content = f.read()
        st.components.v1.html(html_content, height=750, width=1000)
    
    # Node/Edge modification section
    st.subheader("Add Nodes and Edges")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("Add a new node")
        new_node = st.text_input("Node name:")
        if st.button("Add Node", key="add_node"):
            if new_node and new_node not in st.session_state.added_nodes:
                st.session_state.added_nodes.append(new_node)
                st.success(f"Added node: {new_node}")
                # Regenerate the network
                if st.session_state.df is not None:
                    net, _ = draw_network(st.session_state.df)
                    net.save_graph("network.html")
                    with open("network.html", 'r', encoding='utf-8') as f:
                        html_content = f.read()
                    st.rerun()
    
    with col2:
        st.write("Add a new edge")
        
        # Get all available nodes
        all_nodes = st.session_state.network_nodes + st.session_state.added_nodes if 'network_nodes' in st.session_state else st.session_state.added_nodes
        
        if all_nodes:
            source = st.selectbox("Source node:", all_nodes, key="source")
            target = st.selectbox("Target node:", [n for n in all_nodes if n != source], key="target")
            edge_type = st.radio("Edge type:", ["Activation (1)", "Inhibition (2)"], key="edge_type")
            
            edge_type_val = 1 if edge_type == "Activation (1)" else 2
            
            if st.button("Add Edge", key="add_edge"):
                # Check if an edge between these nodes already exists
                existing_edge_index = None
                for i, (s, t, _) in enumerate(st.session_state.added_edges):
                    if s == source and t == target:
                        existing_edge_index = i
                        break
                
                # If edge exists, replace it
                if existing_edge_index is not None:
                    st.session_state.added_edges[existing_edge_index] = (source, target, edge_type_val)
                    st.success(f"Updated edge: {source} -> {target} to {edge_type}")
                # Otherwise add as new edge
                else:
                    st.session_state.added_edges.append((source, target, edge_type_val))
                    st.success(f"Added edge: {source} -> {target} ({edge_type})")
                
                # Regenerate the network
                if st.session_state.df is not None:
                    net, _ = draw_network(st.session_state.df)
                    net.save_graph("network.html")
                    with open("network.html", 'r', encoding='utf-8') as f:
                        html_content = f.read()
                    st.rerun()
        else:
            st.warning("No nodes available. Please add nodes first.")
    
# Display current added nodes and edges
if st.session_state.added_nodes or st.session_state.added_edges:
    st.subheader("Added Elements")
    
    # Show added nodes
    if st.session_state.added_nodes:
        st.write("Added nodes:", ", ".join(st.session_state.added_nodes))
    
    # Show added edges as a DataFrame
    if st.session_state.added_edges:
        st.write("Added edges:")
        
        # Create a DataFrame to display edges in tabular format
        edge_data = []
        for source, target, edge_type in st.session_state.added_edges:
            edge_data.append({
                "Source": source,
                "Target": target,
                "Type": edge_type,
            })
        
        if edge_data:
            edge_df = pd.DataFrame(edge_data)
            st.dataframe(edge_df[["Source", "Target", "Type"]])
            
            # Add download button for the edges
            csv = edge_df.to_csv(sep='\t', index=False).encode('utf-8')
            st.download_button(
                "Download added edges as TSV",
                csv,
                "added_edges.tsv",
                "text/tab-separated-values",
                key="download-edges"
            )
    
    # Button to create TSV file with all nodes and edges (both original and added)
    if st.button("Generate Complete TSV"):
        # Start with original data
        complete_data = []
        
        if st.session_state.df is not None:
            # Add original edges from the dataframe
            for _, row in st.session_state.df.iterrows():
                complete_data.append({
                    "Source": row["Source"],
                    "Target": row["Target"],
                    "Type": row["Type"],
                    "Citation": row.get("Citation", "") if "Citation" in row else ""
                })
        
        # Add newly added edges
        for source, target, edge_type in st.session_state.added_edges:
            complete_data.append({
                "Source": source,
                "Target": target,
                "Type": edge_type,
                "Citation": ""
            })
        
        # Create DataFrame and show
        if complete_data:
            complete_df = pd.DataFrame(complete_data)
            st.write("Complete network data:")
            st.dataframe(complete_df)
            
            # Add download button for the complete dataset
            csv = complete_df.to_csv(sep='\t', index=False).encode('utf-8')
            st.download_button(
                "Download complete network as TSV",
                csv,
                "complete_network.tsv",
                "text/tab-separated-values",
                key="download-complete"
            )

# Instructions
st.markdown("""
### Instructions:
1. Choose one of the following:
   - Upload a TSV file
   - Paste data directly
   - Use predefined datasets (TRRUST or KEGG Human)
2. To find specific gene interactions:
   - Go to the "Predefined Datasets" tab
   - Enter source and target gene names
   - Click "Find Interaction"
3. Click "Generate Network" to visualize
4. Add new nodes and edges as needed
5. Network will update when you regenerate it

### About the predefined datasets:
- **TRRUST Human**: Transcriptional regulatory relationships in humans
- **KEGG Human**: Gene regulatory relationships from KEGG database
""")
