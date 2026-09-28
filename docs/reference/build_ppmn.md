# Build the foundational Posterior Predictive Meta-analytic Network (PPMN)

Build the foundational Posterior Predictive Meta-analytic Network (PPMN)

## Usage

``` r
build_ppmn(
  formula,
  data,
  txt_size = 3,
  vertex_size = 15,
  vertex_width = 0.21,
  vertex_height = 0.5,
  fun_width = 0.5,
  fun_height = 0.5,
  fun_lwd = 0.4,
  vertex_lwd = 0.4,
  arrow_size = 2,
  arrow_offset = 5.5,
  layout_type = "auto",
  circular = F
)
```

## Arguments

- formula:

  An list of expressions that defines the network structure.

- data:

  A list with estimates and standard errors to be synthesized.

- txt_size:

  Size of the text with a visualization of the graph.

- vertex_size:

  Diameter of the vertices in the DAG.

- vertex_width:

  Width of the variable vertices in the Bipartite graph.

- vertex_height:

  Height of the variable vertices in the Bipartite graph.

- fun_width:

  Width of the vertices representing the edge functions in the Bipartite
  graph.

- fun_height:

  Height of the vertices representing the edge functions in the
  Bipartite graph.

- fun_lwd:

  Width of the edge function vertices.

- vertex_lwd:

  Width of the outlines variable vertices.

- arrow_size:

  Size of the arrow lines.

- arrow_offset:

  Distance of the arrows from the vertices.

- layout_type:

  Type of layout (see ?igraph::layout).

- circular:

  Circular layout or not (default circular = FALSE).
