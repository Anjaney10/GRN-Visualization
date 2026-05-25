# GRN-Visualization

GRN-Visualization is an interactive Python application for exploring and visualizing gene regulatory networks (GRNs). It is designed for researchers and enthusiasts working with regulatory genomics, providing easy access to curated biological interaction databases and flexible tools to add, remove, and analyze network components.

## Features

- Visualize gene regulatory networks interactively.
- Query known regulatory links from integrated public databases.
- Add or remove nodes and edges dynamically within the network graph.
- Compare and cross-reference data across multiple databases.
- Intuitive user interface for network exploration and editing.

## Databases Used

This tool integrates regulatory interactions from the following public resources:

- **[RegNetwork](https://doi.org/10.1093/database/bav095)**  
  RegNetwork is a comprehensive repository that collates regulatory relationships among transcription factors (TFs), microRNAs (miRNAs), and target genes, extracted from curated publications and expert reviews.

- **[TRRUST](https://doi.org/10.1038/srep11432)**  
  TRRUST is an expert-curated database providing human and mouse transcriptional regulatory networks, focusing on TF-target regulatory relationships verified by literature evidence.

> *Planned*: Additional and more updated databases are being considered for future releases to enhance cross-validation and analytical power.

## Visualization Methods

The core of GRN-Visualization leverages popular Python libraries to render and interact with GRNs:

- **NetworkX**: Used for constructing and managing the underlying graph data structure representing the regulatory network.
- **Matplotlib**: For static plotting of the network graphs.

Typical steps in the visualization pipeline include:

1. **Data Acquisition:** Import gene regulatory data from RegNetwork and TRRUST based on the user's query.
2. **Graph Construction:** Nodes represent genes or regulatory elements; edges represent regulatory relationships.
3. **User Interaction:** Users can add/delete nodes and edges via the application interface.
4. **Rendering:** The graph is visualized with options for layout customization, highlighting, and filtering based on user input or node/edge attributes.

## Installation

```bash
git clone https://github.com/Anjaney10/GRN-Visualization.git
cd GRN-Visualization
pip install -r requirements.txt
python app.py
```

## Usage

- Launch the application as described above.
- Use the interface to search, visualize, and modify GRNs.
- Query regulatory relationships by selecting genes or transcription factors of interest.
- Export results or images of network graphs as needed.

## Contributing

We welcome contributions! Please open issues for bug reports, suggestions, or submit pull requests with improvements.

## License

[MIT License](LICENSE)

---

**Contact:** For questions, reach out via [GitHub Issues](https://github.com/Anjaney10/GRN-Visualization/issues).

````
