documentHTML <- "<div>
    <h1>
        Documentation
    </h1>
    <h2>
        Procedures
    </h2>
    <p>
        Step 1: Upload your gene list and Gene Matrix Transposed file (.gmt) file in the Input Files panel
    </p>
    <p>
        Step 2: Click “Run ActivePathways” button and then specify parameters (Cutoff, Significant, Merge Method, Correction Method, Geneset Filter, Overlap measures) in the pop-up dialogue. Click OK to proceed.
    </p>
    <p>
        Step3: A enrichment map visualization of result will be displayed on the right panel.
    </p>
    <h2>
        Features
    </h2>
    <h3>
        Example:
    </h3>
    <p>
Click “Run Example” button and then a small subset of example enrichment map visualization from data in ActivePathways package is displayed.
</p>
<h3>
Node & edge manipulation
</h3>
<p>
You can select the node and position it freely around the canvas. You can zoom-in by scrolling your mouse. Our nodes also support multiple nodes selection. By pressing command (on mac os) or ctrl (on windows), users are allowed to select multiple nodes and dragging them all together in one piece
</p>
<h3>
Adjustment
</h3>
<p>
You can click the “Adjust” button and move to “Adjust” panel.
</p>
<h3>
Font size
</h3>
<p>
You can move the slider to adjust the font size of node labels
</p>
<h3>
Edge cutoff (similarity) & value
</h3>
<p>
You can filter edges based on value selected. By selecting the edge cutoff value, it will remove edges with overlap sore less than cutoff value. The edge value is calculated based on overlap score ranges from 0 to 1. The default value is 0.25, meaning it will remove edges with overlap score less than 0.25. You can also input a specific cutoff value in the input field.
</p>
<h3>
Highlight Neighbors
</h3>
<p>
You can check the box of highlight neighbors. In doing so, when you click the nodes in the network, only the node itself and its direct neighbors are highlighted.
</p>
<h3>
Navigation buttons
</h3>
<p>
In the bottom of the network, you can find several navigation buttons, allowing you to have more precise ways for movement, zoom-in and zoom-out. You can turn off this feature by select “hide navigation buttons” checkbox.
</p>
<h3>
Add node & edge
</h3>
<p>
You can click the “enable manipulation” to add node or edge to the graph.
</p>
<h3>
Gene set information
</h3>
<p>
You can click the node and an information panel will show up below. In the information panel, the id, name, and genes inside this gene set will be displayed.
</p>
<h1>
</div>"