package ds.finalProject;

import java.util.Arrays;

public class Main {
    public static void main(String[] args) {
        Graph<String> graph = new Graph<>();
        graph.loadGraphFromFile("Graph.txt");

        graph.BFS("M");
    }
}