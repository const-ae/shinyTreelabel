import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";

export class D3TreeSelector {
  // Declare the chart dimensions and margins.
  constructor(id, width = 1400, marginTop = 20, marginRight = 20, marginBottom = 30, marginLeft = 40,
              dx = 15){
    this.id = id
    this.width = width;
    this.marginTop = marginTop;
    this.marginRight = marginRight;
    this.marginBottom = marginBottom;
    this.marginLeft = marginLeft;
    this.dx = dx;

    this.svg = null;
    this.tree = null;
    this.gNode = null;
    this.gLink = null;
    this.root = null;
    this.collapsed_elements = new Set()


    Shiny.addCustomMessageHandler(this.id + "-firstTreeFullData", (message) => {
      this.root = d3.hierarchy(message);

      // Rows are separated by dx pixels, columns by dy pixels. These names can be counter-intuitive
      // (dx is a height, and dy a width). This because the tree must be viewed with the root at the
      // “bottom”, in the data domain. The width of a column is based on the tree’s height.
      const dy = (width - marginRight - marginLeft) / (1 + this.root.height);

      // Define the tree layout and the shape for links.
      this.tree = d3.tree().nodeSize([dx, dy]);

      // Create the SVG container, a layer for the links and a layer for the nodes.
      this.svg = d3.create("svg")
          .attr("width", width)
          .attr("height", dx)
          .attr("viewBox", [-marginLeft, -marginTop, width, dx])
          .attr("style", "max-width: 100%; height: auto; font: 12px sans-serif; user-select: none;");

      this.gLink = this.svg.append("g")
          .attr("fill", "none")
          .attr("stroke", "#555")
          .attr("stroke-opacity", 0.4)
          .attr("stroke-width", 1.5);

      this.gNode = this.svg.append("g")
          .attr("cursor", "pointer")
          .attr("pointer-events", "all");

      this.root.descendants().forEach((d, i) => {
        d.id = d.data.name;
        d._children = d.children;
        if(this.collapsed_elements.has(d.id)){
          d.children = null
        }
      });

      this.updateTree();

      // Append the SVG element.
      console.log("Appending svg element to: " + this.id + "-d3tree_holder");
      document.getElementById(this.id + "-d3tree_holder").append(this.svg.node());
    });


    Shiny.addCustomMessageHandler(this.id + "-treeFullData", (message) => {
     console.log("Update Tree")
     console.log(message);
     this.root = d3.hierarchy(message);

     this.root.descendants().forEach((d, i) => {
        d.id = d.data.name;
        d._children = d.children;
        if(this.collapsed_elements.has(d.id)){
          d.children = null
        }
      });

     this.updateTree();
    });
  }

  treeDimension(root){
    let left = root;
    let right = root;
    this.root.eachBefore(node => {
      if (node.x < left.x) left = node;
      if (node.x > right.x) right = node;
    });

    const height = right.x - left.x + this.marginTop + this.marginBottom;
    return {"start": left.x, "height": height};
  }

  updateTree(collapsing = false) {
    const nodes = this.root.descendants().reverse();
    const links = this.root.links();
    const duration = 250;
    const diagonal = d3.linkHorizontal().x(d => d.y).y(d => d.x);

    // Compute the new tree layout.
    this.tree(this.root);
    const tree_dim = this.treeDimension(this.root);

    // Transition of the svg container
    const transition = this.svg
      .transition()
      .duration(duration)
      .attr("height", tree_dim.height)
      .attr("viewBox", [-this.marginLeft, tree_dim.start - this.marginTop, this.width, tree_dim.height])
      .tween("resize", window.ResizeObserver ? null : () => () => this.svg.dispatch("toggle"));

    const node = this.gNode.selectAll("g")
      .data(nodes, d => d.id)
      .join(
        (enter) => {
          let group = enter.append('g')
            .style("fill-opacity", 0)
            .style("stroke-opacity", 0)
            .attr("transform", d => `translate(2000,${d.x})`)
            .on("click", (event, d) => {
               if(d.data.selectable){
                  if(event?.altKey){
                    if(d.children){
                      this.collapsed_elements.add(d.id);
                      d.children = null
                    }else{
                      this.collapsed_elements.delete(d.id);
                      d.children = d._children
                    }
                    this.updateTree(collapsing = true);
                  }else{
                    Shiny.setInputValue(this.id + "-d3TreeClick", d.data.name, {priority: "event"});
                  }
               }
            })
            .on("mouseover",  (event, d) => {
              if(d.data.selectable){
                  d3.select(event.currentTarget)
                    .select('circle')
                    .attr("r", 10);
              }
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
      .style("fill-opacity", d => {
        return d.data.selectable ? 1 : 0.2
      })
      .style("stroke-opacity", d => {
        return d.data.selectable ? 1 : 0.2
      });

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




    const link = this.gLink.selectAll("path")
      .data(links, d => d.source.id + "-->" + d.target.id)
      .join(
        (enter) => {
          return enter.append("path")
            .style("fill-opacity", 0)
            .style("stroke-opacity", 0);
            // .attr("d", d => {
        },
        (update) => {
          if(collapsing){
            return update
            .style("fill-opacity", 0)
            .style("stroke-opacity", 0);
          }else{
            return update
          }
        },
        (exit) => {
          return exit
            .remove()
            .style("fill-opacity", 0)
            .style("stroke-opacity", 0);
        }
      )
      .attr("d", diagonal)
      .transition()
      .duration(100)
      .delay(duration)   // Delay until global animation is done.
      .style("stroke-opacity", d => {
        return d.target.data.selectable ? 0.4 : 0.03
      });
  }
}

