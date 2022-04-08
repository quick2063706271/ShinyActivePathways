# ShinyActivePathways
A shiny web server that performs ActivePathways and display results in Enrichment map

# Imports
```
shiny,
ActivePathways,
visNetwork,
shinyWidgets,
shinyscreenshot,
shinyalert,
dplyr,
spsComps,
igraph,
data.table
```

# Run 
Run "Run App" on the R studio to run the app or deploy it on a server

# Procedures
Step 1: Upload your gene list and Gene Matrix Transposed file (.gmt) file in the Input Files panel

Step 2: Click “Run ActivePathways” button and then specify parameters (Cutoff, Significant, Merge Method, Correction Method, Geneset Filter, Overlap measures) in the pop-up dialogue. Click OK to proceed.

Step3: A enrichment map visualization of result will be displayed on the right panel.


# Features
## Example: 
Click “Run Example” button and then a small subset of example enrichment map visualization from data in ActivePathways package is displayed. 

## Node & edge manipulation
You can select the node and position it freely around the canvas. You can zoom-in by scrolling your mouse. Our nodes also support multiple nodes selection. By pressing command (on mac os) or ctrl (on windows), users are allowed to select multiple nodes and dragging them all together in one piece

## Adjustment
You can click the “Adjust” button and move to “Adjust” panel.  
Font size: You can move the slider to adjust the font size of node labels

## Edge cutoff (similarity) & value
You can filter edges based on value selected. By selecting the edge cutoff value, it will remove edges with overlap sore less than cutoff value. The edge value is calculated based on overlap score ranges from 0 to 1. The default value is 0.25, meaning it will remove edges with overlap score less than 0.25. You can also input a specific cutoff value in the input field.

## Highlight Neighbors
You can check the box of highlight neighbors. In doing so, when you click the nodes in the network, only the node itself and its direct neighbors are highlighted.

## Navigation buttons
In the bottom of the network, you can find several navigation buttons, allowing you to have more precise ways for movement, zoom-in and zoom-out. You can turn off this feature by select “hide navigation buttons” checkbox.

## Add node & edge
You can click the “enable manipulation” to add node or edge to the graph.

## Gene set information
You can click the node and an information panel will show up below. In the information panel, the id, name, and genes inside this gene set will be displayed. 
