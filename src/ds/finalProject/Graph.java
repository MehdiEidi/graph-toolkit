package ds.finalProject;

import java.io.*;
import java.util.*;

/**
 * This class is capable of storing a directed weighted graph structure. It uses hashTables for storing the graph.
 * @param <K> type of the vertices of the graph
 */
@SuppressWarnings("unchecked")
public class Graph<K> implements Iterable<HashTable.Node<K, HashTable<K, Double>>> {
    //Each graph has a hashTable in which the keys are the label's of the nodes and the values are hashTables in which store all the direct adjacent
    //nodes. The second hashTable's keys are labels of the nodes and the values are the weight between this node and it's parent.
    private final HashTable<K, HashTable<K, Double>> adjacencyTable = new HashTable<>();

    /**
     * Checks to see if a node already exists in the graph.
     * @param nodeKey the key of the node to look for
     * @return true, if the node exists in the graph
     */
    public boolean nodeExists(K nodeKey) {
        return adjacencyTable.containsKey(nodeKey);
    }

    /**
     * Returns the weight of the edge between startNode and endNode.
     * @param startNodeKey the key of the startNode
     * @param endNodeKey the key of the endNode
     * @return the weight of the edge between startNode and endNode
     */
    public Double getWeight(K startNodeKey, K endNodeKey) {
        HashTable<K, Double> adjacentTable = adjacencyTable.get(startNodeKey);

        if (adjacentTable == null) {
            throw new IllegalArgumentException("Edge from \"" + startNodeKey.toString() + "\" to \"" + endNodeKey.toString() + "\" doesn't exist.");
        }

        try {
            return adjacentTable.get(endNodeKey);
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException("No edges between \"" + startNodeKey + "\" and \"" + endNodeKey + "\"");
        }
    }

    /**
     * Adds a new node to the graph.
     * @param nodeKey the key of the node to be added.
     */
    public void addNode(K nodeKey) {
        if (nodeExists(nodeKey)) {
            throw new IllegalArgumentException("This node already exists: " + nodeKey.toString());
        }

        adjacencyTable.put(nodeKey, new HashTable<>());
    }

    /**
     * Adds a new edge to the graph
     * @param startNodeKey the key of the startNode
     * @param endNodeKey the key of the endNode
     * @param weight the weight of the edge
     */
    public void addEdge(K startNodeKey, K endNodeKey, Double weight) {
        //If any of the startNode or endNode doesn't exist, first we need to add them to the graph.
        if (!nodeExists(startNodeKey)) {
            this.addNode(startNodeKey);
        }

        if (!nodeExists(endNodeKey)) {
            this.addNode(endNodeKey);
        }

        adjacencyTable.get(startNodeKey).put(endNodeKey, weight);
    }

    /**
     * Adds a new edge to the graph. The weight of this edge will be zero.
     * @param startNodeKey the key of the startNode
     * @param endNodeKey the key of the endNode
     */
    public void addEdge(K startNodeKey, K endNodeKey) {
        addEdge(startNodeKey, endNodeKey, 0.0);
    }

    /**
     * Removes a node from the graph.
     * @param nodeKey the key of the node to be removed.
     * @return neighbors of the removed node's start node.
     */
    public HashTable<K, Double> removeNode(K nodeKey) {
        for (HashTable.Node<K, HashTable<K, Double>> currentNode : this) {
            if (currentNode.getValue().containsKey(nodeKey)) {
                currentNode.getValue().remove(nodeKey);
            }
        }

        return adjacencyTable.remove(nodeKey);
    }

    /**
     * Removes an edge from the graph.
     * @param startNodeKey the key of the startNode
     * @param endNodeKey the key of the endNode
     * @return the weight of that edge.
     */
    public Double removeEdge(K startNodeKey, K endNodeKey) {
        return adjacencyTable.get(startNodeKey).remove(endNodeKey);
    }

    /**
     * @param nodeKey the key of the node
     * @return a set containing all the neighbors of the given node.
     */
    public Set<HashTable.Node<K, HashTable<K, Double>>> getNeighbors(K nodeKey) {
        Set<HashTable.Node<K, HashTable<K, Double>>> set = new HashSet<>();

        for (HashTable.Node<K, Double> currentAdjacent : adjacencyTable.getNode(nodeKey).getValue()) {
            set.add(adjacencyTable.getNode(currentAdjacent.getKey()));
        }

        return set;
    }

    /**
     * @param nodeKey the key of the node
     * @return an array containing all the keys of the neighbors of the given node.
     */
    public String[] getNeighborsKeys(K nodeKey) {
        Set<HashTable.Node<K, HashTable<K, Double>>> neighborsKeys = getNeighbors(nodeKey);

        String[] keys = new String[neighborsKeys.size()];
        int index = 0;

        for (HashTable.Node<K, HashTable<K, Double>> currentNode : neighborsKeys) {
            keys[index] = currentNode.getKey().toString();
            index++;
        }

        return keys;
    }

    /**
     * Constructs a graph from the given file. The format of the file is like this: [startNode][a space][endNode (optional)][a space][weight (Optional)]
     * @param fileName the name of the file
     */
    public void loadGraphFromFile(String fileName) {
        try {
            FileInputStream fileInputStream = new FileInputStream(fileName);
            Scanner scanner = new Scanner(fileInputStream);

            //Scan the file line by line. Split the line into a string array.
            //line[0] will be the start node. line[1] will be the end node. line[2] will be the weight of the edge.
            while (scanner.hasNextLine()) {
                String[] line = scanner.nextLine().trim().split("\\s+");

                //If the weight didn't mentioned, we assume it's zero.
                if (line.length == 2) {
                    this.addEdge((K) line[0],(K) line[1], 0.0);
                } else if (line.length < 2) { //If only one node is present in the line, we add that node to the graph.
                    this.addNode((K) line[0]);
                } else {
                    this.addEdge((K) line[0],(K) line[1], Double.parseDouble(line[2]));
                }
            }

            scanner.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Writes the graph to a file. The format of the file will be like this: [startNode][a space][endNode][a space][weight]
     * @param fileName the name of the file
     */
    public void writeGraphToFile(String fileName) {
        File outputFile = new File(fileName);

        try {
            FileOutputStream fileOutputStream = new FileOutputStream(outputFile);
            BufferedWriter bufferedWriter = new BufferedWriter(new OutputStreamWriter(fileOutputStream));

            try {
                //Iterates the graph.
                for (HashTable.Node<K, HashTable<K, Double>> currentNode : this) {
                    for (HashTable.Node<K, Double> currentAdjacent : currentNode.getValue()) {
                        bufferedWriter.write(currentNode.getKey() + " " + currentAdjacent.getKey() + " " + currentAdjacent.getValue());
                        bufferedWriter.newLine();
                    }
                }

                bufferedWriter.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    //Because a graph may contain some sub graphs which are not connected to each other, so, if we do DFS on one of those sub graphs, others
    //will remain unvisited. So, we need to loop over all nodes and perform DFS on them if they are not already visited. That's the reason
    //we have a helper DFS method here.
    private void DFSMain(K nodeKey, Set<HashTable.Node<K, HashTable<K, Double>>> visitedNodes) {
        Stack<HashTable.Node<K, HashTable<K, Double>>> stack = new Stack<>();

        stack.push(adjacencyTable.getNode(nodeKey));

        while (!stack.isEmpty()) {
            HashTable.Node<K, HashTable<K, Double>> currentNode = stack.pop();

            if (!visitedNodes.contains(currentNode)) {
                visitedNodes.add(currentNode);

                System.out.print(currentNode.getKey() + "  ");

                for (HashTable.Node<K, Double> currentAdjacent : currentNode.getValue()) {
                    stack.push(adjacencyTable.getNode(currentAdjacent.getKey()));
                }
            }
        }
    }

    /**
     * Traverses the graph using DFS algorithm.
     * @param nodeKey the key of the node that we want to start the search with.
     */
    public void DFS(K nodeKey) {
        Set<HashTable.Node<K, HashTable<K, Double>>> visitedNodes = new HashSet<>(); //For storing the nodes that are already visited. We don't want to revisit the nodes that are already visited.
        DFSMain(nodeKey, visitedNodes);

        for (HashTable.Node<K, HashTable<K, Double>> currentNode : this) {
            if (!visitedNodes.contains(currentNode)) {
                DFSMain(currentNode.getKey(), visitedNodes);
            }
        }
    }

    //Because a graph may contain some sub graphs which are not connected to each other, so, if we do BFS on one of those sub graphs, others
    //will remain unvisited. So, we need to loop over all nodes and perform BFS on them if they are not already visited. That's the reason
    //we have a helper BFS method here.
    private void BFSMain(K nodeKey, Set<HashTable.Node<K, HashTable<K, Double>>> visitedNodes) {
        Queue<HashTable.Node<K, HashTable<K, Double>>> queue = new ArrayDeque<>();

        queue.add(adjacencyTable.getNode(nodeKey));

        while (!queue.isEmpty()) {
            HashTable.Node<K, HashTable<K, Double>> currentNode = queue.poll();

            if (!visitedNodes.contains(currentNode)) {
                visitedNodes.add(currentNode);

                System.out.print(currentNode.getKey() + "  ");

                for (HashTable.Node<K, Double> currentAdjacent : currentNode.getValue()) {
                    queue.add(adjacencyTable.getNode(currentAdjacent.getKey()));
                }
            }
        }
    }

    /**
     * Traverses the graph using BFS algorithm.
     * @param nodeKey the key of the node that we want to start the search with.
     */
    public void BFS(K nodeKey) {
        Set<HashTable.Node<K, HashTable<K, Double>>> visitedNodes = new HashSet<>(); //For storing the nodes that are already visited. We don't want to revisit the nodes that are already visited.
        BFSMain(nodeKey, visitedNodes);

        for (HashTable.Node<K, HashTable<K, Double>> currentNode : this) {
            if (!visitedNodes.contains(currentNode)) {
                BFSMain(currentNode.getKey(), visitedNodes);
            }
        }
    }

    //Helper method for topological sort algorithm. It repeats the process of topological sort on all the neighbors of the given node.
    private void topologicalSortHelper(Set<HashTable.Node<K, HashTable<K, Double>>> visitedNodes, Stack<HashTable.Node<K, HashTable<K, Double>>> topologicalOrder, HashTable.Node<K, HashTable<K, Double>> currentNode) {
        visitedNodes.add(currentNode);

        for (HashTable.Node<K, Double> currentAdjacent : currentNode.getValue()) {
            topologicalSortHelper(visitedNodes, topologicalOrder, adjacencyTable.getNode(currentAdjacent.getKey()));
        }

        //The given node may be already added to the stack. So, we don't want to re add that node to the stack.
        if (!topologicalOrder.contains(currentNode)) {
            topologicalOrder.push(currentNode);
        }
    }

    /**
     * Sorts the graph in the topological order. The order will be stored in a stack.
     * @param topologicalOrder for storing the topological order.
     */
    public void topologicalSort(Stack<HashTable.Node<K, HashTable<K, Double>>> topologicalOrder) {
        Set<HashTable.Node<K, HashTable<K, Double>>> visitedNodes = new HashSet<>(); //For storing the nodes that are already visited. We don't want to revisit the nodes which are already visited.

        //For every node of the graph, if it's not visited already, we do the process of topological sort on it.
        for (HashTable.Node<K, HashTable<K, Double>> current : this) {
            if (!visitedNodes.contains(current)) {
                topologicalSortHelper(visitedNodes, topologicalOrder, current);
            }
        }
    }

    /**
     * Sorts the graph in the topological, stores the order in a stack and prints the stack.
     */
    public void printTopologicalSort() {
        Stack<HashTable.Node<K, HashTable<K, Double>>> topologicalOrder = new Stack<>(); //Stores the topological order of the graph.

        topologicalSort(topologicalOrder);

        while (!topologicalOrder.isEmpty()) {
            System.out.print(topologicalOrder.pop().getKey() + "   ");
        }
    }

    /**
     * Finds shortest-path distances from the given source node to all other nodes in the graph using topological sorting. But this algorithm
     * will only work for the directed acyclic graphs. If a graph is not acyclic you can use the dijkstra's algorithm available in the class.
     * @param sourceNodeKey the key of the source node.
     */
    public void topologicalShortestPath(K sourceNodeKey) {
        Stack<HashTable.Node<K, HashTable<K, Double>>> topologicalOrder = new Stack<>(); //For storing the topological order of the graph.
        topologicalSort(topologicalOrder);

        Map<HashTable.Node<K, HashTable<K, Double>>, Double> distances = new HashMap<>(); //For storing entries of a node and it's distance to the source node.

        //Marking all nodes distances to the source node with INFINITY except the source node itself. We must mark the source node's distance to itself with zero.
        for (HashTable.Node<K, HashTable<K, Double>> currentNode : topologicalOrder) {
            if (currentNode.getKey().equals(sourceNodeKey)) {
                distances.put(currentNode, 0.0);
            } else {
                distances.put(currentNode, Double.POSITIVE_INFINITY);
            }
        }

        //Processing the nodes using the topological order stored in the stack.
        while (!topologicalOrder.isEmpty()) {
            HashTable.Node<K, HashTable<K, Double>> currentNode = topologicalOrder.pop();

            if (distances.get(currentNode) != Double.POSITIVE_INFINITY) {
                for (HashTable.Node<K, Double> currentAdjacent : currentNode.getValue()) {
                    if (distances.get(adjacencyTable.getNode(currentAdjacent.getKey())) > (distances.get(currentNode) + currentAdjacent.getValue())) {
                        distances.replace(adjacencyTable.getNode(currentAdjacent.getKey()), distances.get(currentNode) + currentAdjacent.getValue());
                    }
                }
            }
        }

        //Printing the entries of the "distances" map.
        printMap(distances, sourceNodeKey);
    }

    //Utility method for iterating and printing a map. This method is used in the shortest-path distance methods.
    private void printMap(Map<HashTable.Node<K, HashTable<K, Double>>, Double> map, K sourceNodeKey) {
        System.out.println("Shortest-path distances from the node \"" + sourceNodeKey + "\" to all other nodes are listed below:");

        for (Map.Entry<HashTable.Node<K, HashTable<K, Double>>, Double> entry : map.entrySet()) {
           if (!entry.getKey().getKey().equals(sourceNodeKey) && (entry.getValue() != Double.POSITIVE_INFINITY)) {
                System.out.println("\"" + sourceNodeKey + "\"" + " -> \"" + entry.getKey().getKey() + "\" : " + entry.getValue());
            } else if (!entry.getKey().getKey().equals(sourceNodeKey)) {
                System.out.println("We cannot reach from the node \"" + sourceNodeKey + "\" to the node \"" + entry.getKey().getKey() + "\"");
            }
        }
    }

    /**
     * Finds the shortest-path distances from the given source node to all other nodes in the graph using Dijkstra's algorithm.
     * @param sourceNodeKey the key of the source node
     */
    public Map<HashTable.Node<K, HashTable<K, Double>>, Double> shortestPaths(K sourceNodeKey) {
        Map<HashTable.Node<K, HashTable<K, Double>>, Double> distances = new HashMap<>(); //For storing entries of (node, distance).
        Set<HashTable.Node<K, HashTable<K, Double>>> visitedNodes = new HashSet<>(); //Set of the visited nodes. We don't want to revisit a node which is already visited.
        Queue<HashTable.Node<K, HashTable<K, Double>>> queue = new ArrayDeque<>(); //For using the BFS structure, we need a queue to store the adjacent nodes of a node.

        //Marking all nodes distances from the source node with INFINITY except the source node itself. We must mark the source node's distance to itself with zero.
        for (HashTable.Node<K, HashTable<K, Double>> currentNode : this) {
            if (currentNode.getKey().equals(sourceNodeKey)) {
                distances.put(currentNode, 0.0);
            } else {
                distances.put(currentNode, Double.POSITIVE_INFINITY);
            }
        }

        //Using BFS algorithm to store and process the nodes layer by layer.
        queue.add(adjacencyTable.getNode(sourceNodeKey));
        while (!queue.isEmpty()) {
            HashTable.Node<K, HashTable<K, Double>> currentNode = queue.poll();

            if (!visitedNodes.contains(currentNode)) {
                for (HashTable.Node<K, Double> currentAdjacent : currentNode.getValue()) {
                    //We want to update the distances only if the sum of the current node's distance and the weight of that edge is less than the
                    //existing distance stored in the map.
                    if (distances.get(currentNode) + currentAdjacent.getValue() < distances.get(adjacencyTable.getNode(currentAdjacent.getKey()))) {
                        distances.put(adjacencyTable.getNode(currentAdjacent.getKey()), distances.get(currentNode) + currentAdjacent.getValue());
                    }
                }

                //Adding all the neighbors of the current node to the queue to be processed later.
                queue.addAll(this.getNeighbors(currentNode.getKey()));

                visitedNodes.add(currentNode);
            }
        }

        //Printing the entries of the "distances" map.
        printMap(distances, sourceNodeKey);

        return distances;
    }

    /**
     * Prints the shortest-path distance from sourceNode to destNode
     * @param sourceNodeKey the key of the source node
     * @param destNodeKey the key of the destination node
     */
    public void shortestPathBetween(K sourceNodeKey, K destNodeKey) {
        Map<HashTable.Node<K, HashTable<K, Double>>, Double> distances = shortestPaths(sourceNodeKey); //Storing all shortest-path distances.

        System.out.print("Shortest-path distance from \"" + sourceNodeKey+ "\" to " + "\"" + destNodeKey + "\" is: ");

        for (HashTable.Node<K, HashTable<K, Double>> entry : distances.keySet()) {
            if (entry.getKey().equals(destNodeKey)) {
                System.out.print(distances.get(entry));
            }
        }
    }

    /**
     * Stores all parents of a node in a set and returns that set.
     * @param nodeKey key of the node that we are looking for it's parents.
     * @return a set containing all the parents of the given node.
     */
    public Set<HashTable.Node<K, HashTable<K, Double>>> getParents(K nodeKey) {
        Set<HashTable.Node<K, HashTable<K, Double>>> parents = new HashSet<>();

        //Iterating over the graph's nodes and their children.
        for (HashTable.Node<K, HashTable<K, Double>> currentNode : this) {
            for (HashTable.Node<K, Double> currentAdjacent : currentNode.getValue()) {
                if (currentAdjacent.getKey().equals(nodeKey)) {
                    parents.add(currentNode);
                }
            }
        }

        return parents;
    }

    //Get's a map and a color and checks to see if that map contains that color.
    private boolean containsColor(int color, Map<HashTable.Node<K, HashTable<K, Double>>, Integer> nodesColors) {
        for (HashTable.Node<K, HashTable<K, Double>> entry : nodesColors.keySet()) {
            if (nodesColors.get(entry).equals(color)) {
                return true;
            }
        }

        return false;
    }

    /**
     * Colors the graph's nodes in a way that no adjacent node has the same color.
     */
    public void graphColoring() {
        ArrayList<Integer> colors = new ArrayList<>(); //Stores the available colors. We use numbers as colors. For example 3 is a kind of color.
        colors.add(0);

        Map<HashTable.Node<K, HashTable<K, Double>>, Integer> nodesColors = new HashMap<>(); //For storing the (node, color) relations.

        //Iterating over all the nodes in the graph and give a proper color to that node.
        for (HashTable.Node<K, HashTable<K, Double>> currentNode : this) {
            int colorIndex = 0; //Index of the node's color among the colors array.
            nodesColors.put(currentNode, colors.get(0)); //Initially, we give the node color 0.

            Map<HashTable.Node<K, HashTable<K, Double>>, Integer> adjacentColors = new HashMap<>(); //Stores the colors of the adjacent nodes of the current node.

            for (HashTable.Node<K, Double> currentAdjacent : currentNode.getValue()) {
                //We should first check to see the current adjacent color is given a color or not (we don't care about colorless nodes).
                if (nodesColors.containsKey(adjacencyTable.getNode(currentAdjacent.getKey()))) {
                    adjacentColors.put(adjacencyTable.getNode(currentAdjacent.getKey()), nodesColors.get(adjacencyTable.getNode(currentAdjacent.getKey())));

                    //If the both nodes have the same color.
                    if (nodesColors.get(adjacencyTable.getNode(currentAdjacent.getKey())).equals(nodesColors.get(currentNode))) {
                        while (nodesColors.get(adjacencyTable.getNode(currentAdjacent.getKey())).equals(nodesColors.get(currentNode))) {
                            //If we are out of available colors, we should introduce a new color.
                            if (colorIndex >= colors.size() - 1) {
                                colors.add(colorIndex + 1);
                            }
                            nodesColors.replace(currentNode, colors.get(++colorIndex));
                        }
                    }
                }
            }

            colorIndex = 0;

            Set<HashTable.Node<K, HashTable<K, Double>>> parents = getParents(currentNode.getKey()); //Stores all the parents of the currentNode.

            //Iterating over all parents of the currentNode and check to see if it has the same color with the currentNode. If so, we should either pick another color or introduce a new color if we are out of available colors.
            for (HashTable.Node<K, HashTable<K, Double>> currentParent : parents) {
                if (nodesColors.containsKey(currentParent) && nodesColors.get(currentNode).equals(nodesColors.get(currentParent))) {
                    //While the parent's color is same as currentNode's color or the color that we picked now is the same as one of the adjacent node's color, we should pick a new color.
                    while (nodesColors.get(currentParent).equals(nodesColors.get(currentNode)) || containsColor(colors.get(colorIndex), adjacentColors)) {
                        if (colorIndex >= colors.size() - 1) {
                            colors.add(colorIndex + 1);
                        }
                        nodesColors.replace(currentNode, colors.get(++colorIndex));
                    }
                }
            }
        }

        //Printing the (node, color) entries.
        for (HashTable.Node<K, HashTable<K, Double>> entry : nodesColors.keySet()) {
            System.out.println("\"" + entry.getKey() + "\" color: " + nodesColors.get(entry));
        }
    }

    //Fills the stack using a DFS on the graph starting from the given node.
    private void fillOrder(HashTable.Node<K, HashTable<K, Double>> currentNode, Set<HashTable.Node<K, HashTable<K, Double>>> visitedNodes, Stack<HashTable.Node<K, HashTable<K, Double>>> sequenceOfNodes) {
        visitedNodes.add(currentNode);

        for (HashTable.Node<K, Double> currentAdjacent : currentNode.getValue()) {
            if (!visitedNodes.contains(adjacencyTable.getNode(currentAdjacent.getKey()))) {
                fillOrder(adjacencyTable.getNode(currentAdjacent.getKey()), visitedNodes, sequenceOfNodes);
            }
        }

        sequenceOfNodes.push(currentNode);
    }

    /**
     * Reverses the graph's edge's directions. An edge from "a" to "b" will become an edge from "b" to "a"
     * @return the reversed (transposed) graph.
     */
    public Graph<K> transpose() {
        Graph<K> transposedGraph = new Graph<>();

        for (HashTable.Node<K, HashTable<K, Double>> currentNode : this) {
            for (HashTable.Node<K, Double> currentAdjacent : currentNode.getValue()) {
                transposedGraph.addEdge(currentAdjacent.getKey(), currentNode.getKey(), currentAdjacent.getValue());

            }
        }

        return transposedGraph;
    }

    /**
     * Prints the strongly connected components of the graph using Kosaraju's algorithm.
     */
    public void SCC() {
        Stack<HashTable.Node<K, HashTable<K, Double>>> sequenceOfNodes = new Stack<>(); //For storing sequence of all nodes using a DFS.
        Set<HashTable.Node<K, HashTable<K, Double>>> visitedNodes = new HashSet<>(); //For storing the nodes that are already visited and processed.

        //Performing DFS on all the unvisited nodes of the graph.
        for (HashTable.Node<K, HashTable<K, Double>> currentNode : this) {
            if (!visitedNodes.contains(currentNode)) {
                fillOrder(currentNode, visitedNodes, sequenceOfNodes);
            }
        }

        //Getting the reversed (transposed) graph of the current graph.
        Graph<K> transposedGraph = transpose();

        visitedNodes.clear();

        System.out.println("Strongly connected components of the graph are: (all the keys listed in one line are among one component)");

        //Performing DFS on the nodes of the transposed graph. The order of nodes in which we are performing DFS on, is based on the order of nodes stored in the stack.
        while (!sequenceOfNodes.isEmpty()) {
            HashTable.Node<K, HashTable<K, Double>> currentNode = sequenceOfNodes.pop();

            if (!visitedNodes.contains(currentNode)) {
                transposedGraph.DFSMain(currentNode.getKey(), visitedNodes);
                System.out.println();
            }
        }
    }

    /**
     * @return an iterator of the nodes of the graph.
     */
    public Iterator<HashTable.Node<K, HashTable<K, Double>>> iterator() {
        return new Iterator<>() {
            private HashTable.Node<K, HashTable<K, Double>> currentNode = null;
            private int counter = 0;
            private int tableIndex = -1;

            @Override
            public boolean hasNext() {
                return adjacencyTable.getOccupied() > counter;
            }

            private HashTable.Node<K, HashTable<K, Double>> findNextIndex() {
                tableIndex++;
                for (int i = tableIndex; i < adjacencyTable.getTable().length; i++) {
                    if (adjacencyTable.getTable()[i] == null) {
                        tableIndex++;
                    } else {
                        break;
                    }
                }

                return adjacencyTable.getTable()[tableIndex];
            }

            @Override
            public HashTable.Node<K, HashTable<K, Double>> next() {
                HashTable.Node<K, HashTable<K, Double>> current = currentNode;

                if (currentNode == null) {
                    current = findNextIndex();
                    currentNode = current;
                    counter++;
                    return current;
                } else if (current.getNext() == null) {
                    current = findNextIndex();
                    currentNode = current;
                    counter++;
                    return current;
                } else {
                    current = current.getNext();
                    currentNode = current;
                    counter++;
                    return current;
                }
            }
        };
    }
}