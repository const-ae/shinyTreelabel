import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";

// Declare the chart dimensions and margins.
const width = 1400;
const marginTop = 20;
const marginRight = 20;
const marginBottom = 30;
const marginLeft = 40;

const dx = 15;
const diagonal = d3.linkHorizontal().x(d => d.y).y(d => d.x);
var svg = null;
var tree = null;
var gNode = null;
var gLink = null;

function treeDimension(root){
  let left = root;
  let right = root;
  root.eachBefore(node => {
    if (node.x < left.x) left = node;
    if (node.x > right.x) right = node;
  });

  const height = right.x - left.x + marginTop + marginBottom;
  return {"start": left.x, "height": height};
}

function updateTree(root) {

  const nodes = root.descendants().reverse();
  const links = root.links();
  const duration = 5000;

  // Compute the new tree layout.
  tree(root);
  const tree_dim = treeDimension(root);

  // Transition of the svg container
  const transition = svg
    .transition()
    .duration(duration)
    .attr("height", tree_dim.height)
    .attr("viewBox", [-marginLeft, tree_dim.start - marginTop, width, tree_dim.height])
    .tween("resize", window.ResizeObserver ? null : () => () => svg.dispatch("toggle"));

  const node = gNode.selectAll("g")
    .data(nodes, d => d.id)
    .join(
      (enter) => {
        let group = enter.append('g')
          .style("fill-opacity", 0)
          .style("stroke-opacity", 0)
          .attr("transform", d => `translate(2000,${d.x})`)
          .on("click", (event, d) => {
              if(event?.altKey){
                d.children = d.children ? null : d._children;
                updateTree(root);
              }else{
                Shiny.setInputValue("d3TreeClick", d.data.name);
              }
          })
          .on("mouseover",  (event, d) => {
              d3.select(event.currentTarget)
                .select('circle')
                .attr("r", 10);
          })
          .on("mouseout", (event, d) => {
              d3.select(event.currentTarget)
                .select('circle')
                .attr("r", 5);
          });
        group.append('circle')
          .attr("r", 5)
          // .style('fill', d => d._children ? "#555" : "#999")
        group.append('text')
          .attr("dy", "0.31em")
          .attr("x", d => d._children ? -6 : 6)
          .attr("text-anchor", d => d._children ? "end" : "start")
          .text(d => d.data.name)
          .style("stroke-linejoin", "round")
          .style("stroke-width", 3)
          .style("stroke", "white")
          .style("paint-order", "stroke")
        return group;
      },
      (update) => {
        return update;
      },
      (exit) => {
        return exit
          .transition(transition)
          .remove()
          .style("fill-opacity", 0)
          .style("stroke-opacity", 0);
      }
    )
    .transition(transition)
    .attr("transform", d => `translate(${d.y},${d.x})`)
    .style("fill-opacity", 1)
    .style("stroke-opacity", 1);

  // Update some features of the circle if it is selected.
  node
    .select("circle")
    .style('fill', d => {
      if(d.data.selected){
        return d.data.selectionColor ? d.data.selectionColor : "green";
      }else{
        return d._children ? "#555" : "#999";
      }
    })
    .attr("r", d => d.data.selected ? 7 : 5);




  const link = gLink.selectAll("path")
    .data(links, d => d.source.id + "-->" + d.target.id)
    .join(
      (enter) => {
        console.log("Path enter size: " + enter.size());
        return enter.append("path")
          .style("fill-opacity", 0)
          .style("stroke-opacity", 0);
          // .attr("d", d => {
          //   const o = {x: d.target.x, y: 2000};
          //   const old_source = {x: d.source.x0, y: d.source.y0};
          //   return diagonal({source: old_source, target: o});
          // })
      },
      (update) => {
        console.log("Path update size: " + update.size());
        return update
      },
      (exit) => {
        return exit
          .remove()
          .style("fill-opacity", 0)
          .style("stroke-opacity", 0);
      }
    )
    .attr("d", diagonal)
    .transition(transition)
    .style("fill-opacity", 1)
    .style("stroke-opacity", 1);

}




Shiny.addCustomMessageHandler("firstTreeFullData", function(message) {
  const root = d3.hierarchy(message);

  // Rows are separated by dx pixels, columns by dy pixels. These names can be counter-intuitive
  // (dx is a height, and dy a width). This because the tree must be viewed with the root at the
  // “bottom”, in the data domain. The width of a column is based on the tree’s height.
  const dy = (width - marginRight - marginLeft) / (1 + root.height);

  // Define the tree layout and the shape for links.
  tree = d3.tree().nodeSize([dx, dy]);

  // Create the SVG container, a layer for the links and a layer for the nodes.
  svg = d3.create("svg")
      .attr("width", width)
      .attr("height", dx)
      .attr("viewBox", [-marginLeft, -marginTop, width, dx])
      .attr("style", "max-width: 100%; height: auto; font: 12px sans-serif; user-select: none;");

  gLink = svg.append("g")
      .attr("fill", "none")
      .attr("stroke", "#555")
      .attr("stroke-opacity", 0.4)
      .attr("stroke-width", 1.5);

  gNode = svg.append("g")
      .attr("cursor", "pointer")
      .attr("pointer-events", "all");

  root.descendants().forEach((d, i) => {
    d.id = d.data.name;
    d._children = d.children;
  });

  updateTree(root);

  // Append the SVG element.
  d3tree_holder.append(svg.node());
});


Shiny.addCustomMessageHandler("treeFullData", function(message) {
 console.log("Update Tree")
 console.log(message);
 const root = d3.hierarchy(message);

 root.descendants().forEach((d, i) => {
    d.id = d.data.name;
    d._children = d.children;
  });

 updateTree(root);
});
