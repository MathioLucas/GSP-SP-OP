

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import json
import argparse
import sys
from pathlib import Path

class GraphVisualizer:
    def __init__(self):
        self.graph = nx.Graph()
        self.pos = None
        self.classification = None
        
    def load_graph_from_txt(self, filename):
       
        self.graph.clear()
        
        with open(filename, 'r') as f:
            # First line: vertices and edges
            n_vertices, n_edges = map(int, f.readline().strip().split())
            
            # Add vertices
            self.graph.add_nodes_from(range(n_vertices))
            
            # Add edges
            for _ in range(n_edges):
                u, v = map(int, f.readline().strip().split())
                self.graph.add_edge(u, v)
        
        print(f"Loaded graph: {n_vertices} vertices, {n_edges} edges")
        return self.graph
    
    def load_graph_from_json(self, filename):
       
        with open(filename, 'r') as f:
            data = json.load(f)
        
        self.graph.clear()
        n_vertices = data['vertices']
        self.graph.add_nodes_from(range(n_vertices))
        
        # Load adjacency list
        adj_list = data['adjacency_list']
        for vertex_str, neighbors in adj_list.items():
            vertex = int(vertex_str)
            for neighbor in neighbors:
                if vertex < neighbor:  # Avoid duplicate edges
                    self.graph.add_edge(vertex, neighbor)
        
        # Load classification if available
        if 'classification' in data:
            self.classification = data['classification']
        
        print(f"Loaded graph from JSON: {n_vertices} vertices, {self.graph.number_of_edges()} edges")
        return self.graph
    
    def analyze_graph(self):
        """Compute basic graph statistics"""
        if self.graph.number_of_nodes() == 0:
            return {}
        
        analysis = {
            'vertices': self.graph.number_of_nodes(),
            'edges': self.graph.number_of_edges(),
            'density': nx.density(self.graph),
            'is_connected': nx.is_connected(self.graph),
            'avg_degree': 2 * self.graph.number_of_edges() / self.graph.number_of_nodes() if self.graph.number_of_nodes() > 0 else 0,
            'max_degree': max(dict(self.graph.degree()).values()) if self.graph.number_of_nodes() > 0 else 0,
            'min_degree': min(dict(self.graph.degree()).values()) if self.graph.number_of_nodes() > 0 else 0,
        }
        
        # Additional analysis for connected graphs
        if analysis['is_connected']:
            analysis['diameter'] = nx.diameter(self.graph)
            analysis['avg_shortest_path'] = nx.average_shortest_path_length(self.graph)
            analysis['clustering'] = nx.average_clustering(self.graph)
        
        return analysis
    
    def compute_layout(self, layout_type='spring'):
        """Compute node positions for visualization"""
        if layout_type == 'spring':
            self.pos = nx.spring_layout(self.graph, k=1, iterations=50)
        elif layout_type == 'circular':
            self.pos = nx.circular_layout(self.graph)
        elif layout_type == 'planar':
            if nx.is_planar(self.graph):
                self.pos = nx.planar_layout(self.graph)
            else:
                print("Graph is not planar, using spring layout")
                self.pos = nx.spring_layout(self.graph)
        else:
            self.pos = nx.spring_layout(self.graph)
        
        return self.pos
    
    def visualize_basic(self, title="Graph Visualization", save_path=None):
        """Create basic graph visualization"""
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        if self.pos is None:
            self.compute_layout()
        
        # Node colors based on degree
        node_colors = [self.graph.degree(node) for node in self.graph.nodes()]
        
        # Draw the graph
        nx.draw_networkx_nodes(self.graph, self.pos, 
                              node_color=node_colors, 
                              node_size=300,
                              cmap=plt.cm.viridis,
                              alpha=0.8,
                              ax=ax)
        
        nx.draw_networkx_edges(self.graph, self.pos, 
                              alpha=0.6, 
                              width=1.5,
                              ax=ax)
        
        nx.draw_networkx_labels(self.graph, self.pos, 
                               font_size=8,
                               font_color='white',
                               font_weight='bold',
                               ax=ax)
        
        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.axis('off')
        
        # Add colorbar for degree
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, 
                                  norm=plt.Normalize(vmin=min(node_colors), vmax=max(node_colors)))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('Node Degree', rotation=270, labelpad=15)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Saved visualization to: {save_path}")
        
        return fig, ax
    
    def visualize_with_classification(self, title="Graph Classification", save_path=None):
        """Create visualization highlighting classification results"""
        if not self.classification:
            return self.visualize_basic(title, save_path)
        
        fig, ax = plt.subplots(1, 1, figsize=(14, 10))
        
        if self.pos is None:
            self.compute_layout()
        
        # Color nodes based on classification
        is_gsp = self.classification.get('is_gsp', False)
        is_sp = self.classification.get('is_sp', False)
        is_op = self.classification.get('is_op', False)
        
        # Determine color scheme
        if is_sp and is_op:
            node_color = 'gold'
            class_text = "Both SP and OP"
        elif is_sp and not is_op:
            node_color = 'lightblue'
            class_text = "SP but not OP"
        elif not is_sp and is_op:
            node_color = 'lightgreen'
            class_text = "OP but not SP"
        elif is_gsp:
            node_color = 'lightcoral'
            class_text = "GSP but neither SP nor OP"
        else:
            node_color = 'lightgray'
            class_text = "Not GSP, SP, or OP"
        
        # Draw the graph
        nx.draw_networkx_nodes(self.graph, self.pos, 
                              node_color=node_color, 
                              node_size=400,
                              alpha=0.8,
                              edgecolors='black',
                              linewidths=1,
                              ax=ax)
        
        nx.draw_networkx_edges(self.graph, self.pos, 
                              alpha=0.7, 
                              width=2,
                              ax=ax)
        
        nx.draw_networkx_labels(self.graph, self.pos, 
                               font_size=10,
                               font_weight='bold',
                               ax=ax)
        
        # Add classification info
        info_text = f"Classification: {class_text}\n"
        info_text += f"GSP: {'✓' if is_gsp else '✗'}  "
        info_text += f"SP: {'✓' if is_sp else '✗'}  "
        info_text += f"OP: {'✓' if is_op else '✗'}"
        
        ax.text(0.02, 0.98, info_text, 
                transform=ax.transAxes, 
                fontsize=12,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.axis('off')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Saved classification visualization to: {save_path}")
        
        return fig, ax
    
    def create_statistics_plot(self, analysis, save_path=None):
        """Create a statistics summary plot"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        # Degree distribution
        degrees = [self.graph.degree(node) for node in self.graph.nodes()]
        ax1.hist(degrees, bins=max(1, len(set(degrees))), alpha=0.7, color='skyblue', edgecolor='black')
        ax1.set_xlabel('Degree')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Degree Distribution')
        ax1.grid(True, alpha=0.3)
        
        # Basic statistics
        stats_labels = ['Vertices', 'Edges', 'Avg Degree', 'Max Degree']
        stats_values = [analysis['vertices'], analysis['edges'], 
                       round(analysis['avg_degree'], 2), analysis['max_degree']]
        
        ax2.bar(stats_labels, stats_values, color=['lightcoral', 'lightblue', 'lightgreen', 'gold'])
        ax2.set_ylabel('Count')
        ax2.set_title('Basic Statistics')
        ax2.grid(True, alpha=0.3)
        
        # Connectivity info
        connectivity_labels = ['Density', 'Connected']
        connectivity_values = [analysis['density'], 1 if analysis['is_connected'] else 0]
        
        ax3.bar(connectivity_labels, connectivity_values, color=['orange', 'purple'])
        ax3.set_ylabel('Value')
        ax3.set_title('Connectivity Measures')
        ax3.set_ylim(0, 1.1)
        ax3.grid(True, alpha=0.3)
        
        # Additional metrics (if connected)
        if analysis['is_connected']:
            metrics_labels = ['Diameter', 'Avg Path Length', 'Clustering']
            metrics_values = [analysis.get('diameter', 0), 
                            round(analysis.get('avg_shortest_path', 0), 2),
                            round(analysis.get('clustering', 0), 3)]
            ax4.bar(metrics_labels, metrics_values, color=['red', 'blue', 'green'])
            ax4.set_ylabel('Value')
            ax4.set_title('Graph Metrics')
        else:
            ax4.text(0.5, 0.5, 'Graph is\nDisconnected', 
                    transform=ax4.transAxes, 
                    ha='center', va='center',
                    fontsize=16, fontweight='bold')
            ax4.set_title('Graph Metrics')
        
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Saved statistics plot to: {save_path}")
        
        return fig
    
    def create_comprehensive_report(self, output_dir="visualization_output"):
        """Create a comprehensive visualization report"""
        Path(output_dir).mkdir(exist_ok=True)
        
        analysis = self.analyze_graph()
        
        # Basic visualization
        self.visualize_basic("Basic Graph Structure", 
                           f"{output_dir}/basic_graph.png")
        
        # Classification visualization (if available)
        if self.classification:
            self.visualize_with_classification("Graph Classification", 
                                             f"{output_dir}/classification.png")
        
        # Statistics plot
        self.create_statistics_plot(analysis, 
                                  f"{output_dir}/statistics.png")
        
        # Create different layouts
        layouts = ['spring', 'circular']
        if nx.is_planar(self.graph):
            layouts.append('planar')
        
        for layout in layouts:
            self.compute_layout(layout)
            self.visualize_basic(f"Graph - {layout.title()} Layout", 
                               f"{output_dir}/layout_{layout}.png")
        
        # Save analysis as JSON
        with open(f"{output_dir}/analysis.json", 'w') as f:
            json.dump(analysis, f, indent=2)
        
        print(f"\nComprehensive report saved to: {output_dir}/")
        print("Files created:")
        print("- basic_graph.png")
        if self.classification:
            print("- classification.png")
        print("- statistics.png")
        for layout in layouts:
            print(f"- layout_{layout}.png")
        print("- analysis.json")
        
        return analysis

def main():
    parser = argparse.ArgumentParser(description='Visualize graphs from GSP-SP-OP project')
    parser.add_argument('input_file', help='Input graph file (.txt or .json)')
    parser.add_argument('--output', '-o', default='visualization_output', 
                       help='Output directory for visualization files')
    parser.add_argument('--layout', '-l', default='spring',
                       choices=['spring', 'circular', 'planar'],
                       help='Layout algorithm for graph visualization')
    parser.add_argument('--show', '-s', action='store_true',
                       help='Show plots interactively')
    
    args = parser.parse_args()
    
    visualizer = GraphVisualizer()
    
    # Load graph based on file extension
    if args.input_file.endswith('.json'):
        visualizer.load_graph_from_json(args.input_file)
    else:
        visualizer.load_graph_from_txt(args.input_file)
    
    if visualizer.graph.number_of_nodes() == 0:
        print("Error: No graph loaded or graph is empty")
        sys.exit(1)
    
    # Create comprehensive report
    analysis = visualizer.create_comprehensive_report(args.output)
    
    # Print analysis to console
    print("\n" + "="*50)
    print("GRAPH ANALYSIS SUMMARY")
    print("="*50)
    for key, value in analysis.items():
        if isinstance(value, float):
            print(f"{key.replace('_', ' ').title()}: {value:.4f}")
        else:
            print(f"{key.replace('_', ' ').title()}: {value}")
    
    if args.show:
        plt.show()

if __name__ == "__main__":
    main()
