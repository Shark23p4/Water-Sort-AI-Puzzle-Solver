import java.util.*;

class node {

    String state; // Representing the state of the node
    List<node> children; // List of child nodes
    boolean isGoal; // Indicator for goal state
    String action; // Action that led to this node (for plan generation)
    node parent; // Reference to the parent node for path reconstruction
    int cost; // Cost from parent to node

    // Constructor
    public node(String state) {
        this.state = state;
        this.children = new ArrayList<>();
        this.isGoal = false; // Default to false
        this.action = ""; // Initialize with no action
        this.parent = null; // Initialize with no parent
        this.cost = 0;
    }

    // Method to set the goal status
    public void setGoal(boolean goalStatus) {
        this.isGoal = goalStatus;
    }

    public boolean isGoal() {
        return this.isGoal;
    }
}

class WaterSortSearch extends GenericSearch {

    List<Stack<Character>> bottles;
    int bottleCapacity = 0;
    int numberOfBottles = 0;
    int currcost = 0;

    // Valid colors
    private static final Set<Character> VALID_COLORS = Set.of('r', 'g', 'b', 'y', 'o', 'e');

    // Method to validate initial state format
    public static boolean validateInitialState(String initialState) {
        String[] parts = initialState.split(";");

        // Step 1: Parse the number of bottles and bottle capacity
        int numberOfBottles;
        int bottleCapacity;

        try {
            numberOfBottles = Integer.parseInt(parts[0]); // First part: number of bottles
            bottleCapacity = Integer.parseInt(parts[1]); // Second part: bottle capacity
        } catch (NumberFormatException e) {
            System.out.println("Invalid number of bottles or bottle capacity.");
            return false; // Invalid format for numbers
        }

        // Step 2: Check if the number of bottle definitions matches the number of
        // bottles
        if (parts.length - 2 != numberOfBottles) {
            System.out.println("Mismatch between number of bottles and the actual bottle definitions.");
            return false;
        }

        // Step 3: Validate each bottle's layers
        for (int i = 2; i < parts.length; i++) { // Bottles start from index 2
            String[] layers = parts[i].split(",");

            // Check if the number of layers exceeds the bottle capacity
            if (layers.length > bottleCapacity) {
                System.out.println("Bottle " + (i - 2) + " exceeds bottle capacity.");
                return false;
            }

            // Validate each layer's color and ensure 'e' is only at the end
            boolean encounteredEmpty = false; // Track if we've seen 'e' in this bottle

            for (int j = 0; j < layers.length; j++) {
                char color = layers[j].charAt(0);

                // Validate that the color is one of the allowed values
                if (!VALID_COLORS.contains(color)) {
                    System.out.println("Invalid color '" + color + "' in bottle " + (i - 2));
                    return false;
                }

                // Check that 'e' only appears at the end of the bottle
                if (color == 'e') {
                    encounteredEmpty = true; // We encountered an empty layer
                } else if (encounteredEmpty) {
                    // If a non-'e' layer appears after an 'e', it's an invalid state
                    System.out.println("Invalid state in bottle " + (i - 2) + ": non-empty layer after empty layer.");
                    return false;
                }
            }
        }

        // If all checks passed
        return true;
    }

    private List<Character> stackToList(Stack<Character> stack) {
        return stack.stream().toList();
    }

    // Method to parse the initial state
    public static List<Stack<Character>> parseInitialState(String initialState) {
        String[] parts = initialState.split(";");
        int numberOfBottles = Integer.parseInt(parts[0]);
        List<Stack<Character>> bottles = new ArrayList<>();

        for (int i = 2; i < 2 + numberOfBottles; i++) {
            Stack<Character> bottle = new Stack<>();
            String[] layers = parts[i].split(",");
            for (int j = layers.length - 1; j >= 0; j--) {
                if (!layers[j].equals("e")) {
                    bottle.push(layers[j].charAt(0));
                }
            }
            bottles.add(bottle);
        }
        return bottles;
    }

    // Method to check if a move is valid
    public boolean validMove(List<Stack<Character>> bottles, int i, int j) {
        if (bottles.get(i).isEmpty()) {
            return false; // Bottle i is empty

        }
        if (bottles.get(j).size() == bottleCapacity) {
            return false; // Bottle j is full

        }
        if (bottles.get(j).isEmpty() || bottles.get(i).peek().equals(bottles.get(j).peek())) {
            return true;
        }
        return false;
    }

    // Method to perform a move and return the new state as a string
    public List<Stack<Character>> Pour(int i, int j) {
        List<Stack<Character>> newBottles = deepCopy(bottles);
        char color = newBottles.get(i).peek();
        int layersToMove = 0;

        // Count how many contiguous layers of the same color are on top of bottle i
        while (layersToMove < newBottles.get(i).size()
                && newBottles.get(i).get(newBottles.get(i).size() - 1 - layersToMove) == color) {
            layersToMove++;
        }

        int emptySpaceInJ = bottleCapacity - newBottles.get(j).size();
        currcost = Math.min(layersToMove, emptySpaceInJ);

        // Perform the move by transferring the layers
        for (int k = 0; k < currcost; k++) {
            newBottles.get(j).push(newBottles.get(i).pop());
        }

        return newBottles;
    }

    // Helper method to deep copy the list of bottles (stacks)
    private static List<Stack<Character>> deepCopy(List<Stack<Character>> bottles) {
        List<Stack<Character>> newBottles = new ArrayList<>();
        for (Stack<Character> bottle : bottles) {
            Stack<Character> newBottle = new Stack<>();
            newBottle.addAll(bottle);
            newBottles.add(newBottle);
        }
        return newBottles;
    }

    // Helper method to encode the state into a string
    public String encodeState(List<Stack<Character>> bottles) {
        StringBuilder sb = new StringBuilder();
        for (Stack<Character> bottle : bottles) {
            // Convert the bottle stack to an array to access elements from bottom to top
            Character[] bottleArray = bottle.toArray(new Character[0]);
            int fillLevel = bottleArray.length;

            // Add the colors from bottom to top
            for (int i = 0; i < fillLevel; i++) {
                sb.append(bottleArray[i]);
            }

            // Fill the remaining spaces with "e"s if the bottle isn't full
            for (int i = fillLevel; i < bottleCapacity; i++) {
                sb.append('e');
            }

            sb.append(";"); // Separate bottles
        }
        return sb.toString();
    }

    // Method to parse bottleCapacity from the initial state
    public static int getBottleCapacity(String initialState) {
        String[] parts = initialState.split(";");
        return Integer.parseInt(parts[1]); // Bottle capacity is the second part of the initialState string
    }

    // Check if the current state is a goal state (all bottles contain one color or
    // are empty)
    public static boolean isGoalState(List<Stack<Character>> bottles) {
        for (Stack<Character> bottle : bottles) {
            if (!bottle.isEmpty() && new HashSet<>(bottle).size() > 1) {
                return false;
            }
        }
        return true;
    }

    // Method to generate all possible moves and return a tree of moves
    public node generateMoveTree(String initialState) {
        bottles = parseInitialState(initialState);
        bottleCapacity = getBottleCapacity(initialState); // Derive bottleCapacity from initialState
        numberOfBottles = bottles.size();
        node root = new node(encodeState(bottles));
        Set<String> visitedStates = new HashSet<>();
        generateMoves(root, bottles, visitedStates);
        return root;
    }

    // Recursive method to generate possible moves for a given state, avoiding
    // repeated states and stopping at goal state
    private void generateMoves(node node, List<Stack<Character>> currBottles,
            Set<String> visitedStates) {

        // Check if the current state is a goal state
        if (isGoalState(currBottles)) {
            node.setGoal(true);
            return; // Stop recursion if the goal state is reached
        }

        String currentState = encodeState(currBottles);
        visitedStates.add(currentState); // Mark current state as visited

        // Generate all valid moves and create child nodes
        for (int i = 0; i < numberOfBottles; i++) {
            for (int j = 0; j < numberOfBottles; j++) {
                if (i != j && validMove(currBottles, i, j)) {
                    bottles = currBottles;
                    List<Stack<Character>> newBottles = Pour(i, j);
                    String newState = encodeState(newBottles);

                    // Only proceed if this new state hasn't been visited yet
                    if (!visitedStates.contains(newState)) {
                        node child = new node(newState);
                        child.cost = node.cost + currcost;
                        child.parent = node;
                        child.action = "Pour_" + i + "_" + j;
                        node.children.add(child);

                        // Recurse for further moves
                        generateMoves(child, newBottles, visitedStates);
                    }
                }
            }
        }
    }

    public static String solve(String initialState, String strategy, boolean visualize) {
        validateInitialState(initialState);
        // Parse user input
        WaterSortSearch waterSortSearch = new WaterSortSearch();
        ;
        int heuristicType = 0;
        if (strategy.matches("GR[12]") || strategy.matches("AS[12]")) {
            heuristicType = Character.getNumericValue(strategy.charAt(2));
            strategy = strategy.substring(0, 2); // Trim heuristic type
        }

        long startTime = System.nanoTime(); // Start tracking time
        long startMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory(); // Start tracking
        // memory

        node searchTree = waterSortSearch.generateMoveTree(initialState);
        String solution = "";
        switch (strategy) {
            case "BF":
                solution = BFS(searchTree);
                break;
            case "DF":
                solution = DFS(searchTree);
                break;
            case "ID":
                solution = IDS(searchTree);
                break;
            case "UC":
                solution = UCS(searchTree);
                break;
            case "GR":
                if (heuristicType == 0) {
                    throw new IllegalArgumentException("Invalid strategy: " + strategy);
                }
                solution = GreedySearch(searchTree, heuristicType);
                break;
            case "AS":
                if (heuristicType == 0) {
                    throw new IllegalArgumentException("Invalid strategy: " + strategy);
                }
                solution = AStarSearch(searchTree, heuristicType);
                break;
            default:
                throw new IllegalArgumentException("Invalid strategy: " + strategy);
        }

        long endTime = System.nanoTime(); // End time tracking
        long endMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory(); // End memory tracking

        // Calculate time and memory used
        long timeElapsed = endTime - startTime;
        long memoryUsed = endMemory - startMemory;

        if (visualize) {
            System.out.println("Time elapsed: " + timeElapsed / 1000000 + " ms");
            System.out.println("Memory used: " + memoryUsed / 1024 + " KB");
            String[] solutionParts = solution.split(";");
            System.out.println("plan: " + solutionParts[0]);
            if (solutionParts.length > 1) {
                System.out.println("pathCost: " + solutionParts[1]);
                System.out.println("nodesExpanded: " + solutionParts[2]);
            }
        }

        return solution;
    }

    public static void printTree(node node, int depth) {
        if (node == null) {
            return;
        }

        // Print the current node's state and moves
        System.out.println("Depth " + depth + ": " + node.state);
        for (node child : node.children) {
            System.out.println("    Move: " + child.action + ", Resulting State: " + child.state);
        }

        // Recursively print children
        for (node child : node.children) {
            printTree(child, depth + 1);
        }
    }

    public static void main(String[] args) {
        String initialState = "6;" +
                "4;" +
                "g,g,g,r;" +
                "g,y,r,o;" +
                "o,r,o,y;" +
                "y,o,y,b;" +
                "r,b,b,b;" +
                "e,e,e,e;";

        String searchType = "UC"; // Example input startegy

        boolean visualize = true;

        solve(initialState, searchType, visualize);

    }
}

class GenericSearch {

    private static String generatePlan(node node) {
        // Generate a solution plan based on the node's path from the root
        List<String> moves = new ArrayList<>();
        node current = node;

        // Backtrack to reconstruct the path from the goal to the root
        while (current != null) {
            if (!current.action.isEmpty()) {
                moves.add(current.action); // Collect the action
            }
            current = current.parent; // Move to the parent node
        }

        Collections.reverse(moves); // Reverse to get path from root to goal
        return String.join(",", moves); // Join moves with commas
    }

    private static int calculatePathCost(node node) {

        return node.cost;
    }

    // Helper method to calculate the depth of a node by counting the path to the
    // root
    private static int getNodeDepth(node node) {
        int depth = 0;
        while (node.parent != null) {
            depth++;
            node = node.parent;
        }
        return depth;
    }

    public static String BFS(node searchTree) {
        // BFS Implementation
        int nodesExpanded = 0;
        StringBuilder plan = new StringBuilder();
        int pathCost;

        Queue<node> frontier = new LinkedList<>();
        frontier.add(searchTree);

        while (!frontier.isEmpty()) {
            node currentNode = frontier.poll();
            nodesExpanded++;

            // Check if the current node is the goal
            if (currentNode.isGoal()) {
                plan.append(generatePlan(currentNode)); // Generate plan from goal node
                pathCost = calculatePathCost(currentNode); // Calculate path cost
                return plan.toString() + ";" + pathCost + ";" + nodesExpanded; // Return the result
            }

            // Add all children to the frontier
            for (node neighbor : currentNode.children) {
                neighbor.parent = currentNode; // Set the parent for path reconstruction
                frontier.add(neighbor);
            }
        }

        return "NOSOLUTION"; // Return if no solution is found
    }

    public static String DFS(node searchTree) {
        // DFS Implementation
        int nodesExpanded = 0;
        StringBuilder plan = new StringBuilder();
        int pathCost;

        // Stack for DFS
        Stack<node> frontier = new Stack<>();
        frontier.push(searchTree); // Add the root node to the stack

        while (!frontier.isEmpty()) {
            node currentNode = frontier.pop(); // Pop the node from the stack
            nodesExpanded++; // Track the number of expanded nodes

            // Check if the current node is the goal
            if (currentNode.isGoal()) {
                plan.append(generatePlan(currentNode)); // Generate the plan from goal node
                pathCost = calculatePathCost(currentNode); // Calculate path cost
                return plan.toString() + ";" + pathCost + ";" + nodesExpanded; // Return the solution
            }

            // Create a new stack for the current node's children
            Stack<node> childStack = new Stack<>();
            // Add all children to the new stack
            for (node neighbor : currentNode.children) {
                childStack.push(neighbor); // Add each child node to the child stack
            }

            // Add children to the original stack in reverse order
            while (!childStack.isEmpty()) {
                frontier.push(childStack.pop()); // Pop from childStack and push to frontier (DFS behavior)
            }

        }

        return "NOSOLUTION"; // Return if no solution is found
    }

    public static String IDS(node searchTree) {
        int depth = 0;
        int totalNodesExpanded = 0; // Track total nodes expanded
        String result;

        // Iteratively deepen the depth until the goal is found or the search space is
        // exhausted
        while (true) {
            // Perform depth-limited DFS and return the result including nodes expanded
            result = depthLimitedDFS(searchTree, depth);

            if (result.equals("NOSOLUTION")) {
                return "NOSOLUTION";
            }
            // Parse the result to extract nodesExpanded and other parts
            String[] resultParts = result.split(";");
            int nodesExpandedAtDepth = Integer.parseInt(resultParts[resultParts.length - 1]); // Last part is nodes
                                                                                              // expanded
            totalNodesExpanded += nodesExpandedAtDepth; // Accumulate total nodes expanded

            // If a solution is found or no solution, modify the nodesExpanded part and
            // return it
            if (!result.startsWith("CUTOFF")) {
                // Reconstruct the result string by replacing the nodesExpanded part with
                // totalNodesExpanded
                resultParts[resultParts.length - 1] = String.valueOf(totalNodesExpanded); // Replace nodesExpanded with
                                                                                          // total
                return String.join(";", resultParts); // Rebuild the final string
            }

            // No solution at this depth, increment depth
            depth++;
        }
    }

    // Helper method for depth-limited DFS (DLS)
    private static String depthLimitedDFS(node currentNode, int depthLimit) {
        Stack<node> frontier = new Stack<>();
        frontier.push(currentNode);

        int nodesExpanded = 0; // Track nodes expanded in this iteration
        StringBuilder plan = new StringBuilder();
        boolean cutoffOccurred = false; // Track if a cutoff happens

        while (!frontier.isEmpty()) {
            node node = frontier.pop();
            nodesExpanded++; // Increment nodes expanded for this depth

            // Check if this is the goal node
            if (node.isGoal()) {
                plan.append(generatePlan(node)); // Generate the plan from the goal node
                int pathCost = calculatePathCost(node); // Calculate path cost
                return plan.toString() + ";" + pathCost + ";" + nodesExpanded; // Return solution with nodes expanded
            }

            // Check if the current depth is within the limit
            int nodeDepth = getNodeDepth(node);
            if (nodeDepth < depthLimit) {
                // Add children to the stack if they are not beyond the depth limit
                Stack<node> childStack = new Stack<>();
                for (node neighbor : node.children) {
                    childStack.push(neighbor); // Add each child node to the child stack
                }

                // Add children to the original stack in reverse order
                while (!childStack.isEmpty()) {
                    frontier.push(childStack.pop()); // Pop from childStack and push to frontier (DFS behavior)
                }
            } else if (nodeDepth == depthLimit) {
                // If we hit the depth limit, set the cutoffOccurred flag
                cutoffOccurred = true;
            }
        }

        // If no solution found, return "CUTOFF" and the nodes expanded at this depth
        return cutoffOccurred ? "CUTOFF;" + nodesExpanded : "NOSOLUTION";
    }

    public static String UCS(node searchTree) {
        PriorityQueue<node> frontier = new PriorityQueue<>(Comparator.comparingInt(n -> n.cost));
        Set<String> explored = new HashSet<>();
        frontier.add(searchTree);

        int nodesExpanded = 0;

        while (!frontier.isEmpty()) {
            node currentNode = frontier.poll();
            nodesExpanded++;

            // Check if the current node is the goal
            if (currentNode.isGoal()) {
                String plan = generatePlan(currentNode);
                int pathCost = calculatePathCost(currentNode);
                return plan + ";" + pathCost + ";" + nodesExpanded;
            }

            explored.add(currentNode.state);

            // Expand the node's children
            for (node child : currentNode.children) {
                if (!explored.contains(child.state)) {
                    frontier.add(child);
                }
            }
        }
        return "NOSOLUTION"; // If no solution found
    }

    public static String GreedySearch(node searchTree, int heuristicType) {
        PriorityQueue<node> frontier;
        if (heuristicType == 1) {
            frontier = new PriorityQueue<>(Comparator.comparingInt(n -> n.children.size()));
        } else {
            frontier = new PriorityQueue<>(Comparator.comparing(n -> n.state));
        }
        Set<String> explored = new HashSet<>();
        frontier.add(searchTree);

        int nodesExpanded = 0;

        while (!frontier.isEmpty()) {
            node currentNode = frontier.poll();
            nodesExpanded++;

            // Check if the current node is the goal
            if (currentNode.isGoal()) {
                String plan = generatePlan(currentNode);
                int pathCost = calculatePathCost(currentNode);
                return plan + ";" + pathCost + ";" + nodesExpanded;
            }

            explored.add(currentNode.state);

            // Expand the node's children
            for (node child : currentNode.children) {
                if (!explored.contains(child.state)) {
                    frontier.add(child);
                }
            }
        }
        return "NOSOLUTION"; // If no solution found
    }

    public static String AStarSearch(node searchTree, int heuristicType) {
        PriorityQueue<node> frontier;

        if (heuristicType == 1) {
            // Compare nodes by cost, then by node depth
            frontier = new PriorityQueue<node>(
                    Comparator.comparingInt((node n) -> n.cost)
                            .thenComparingInt((node n) -> n.children.size()));
        } else {
            // Compare nodes by cost, then lexicographically by state
            frontier = new PriorityQueue<node>(
                    Comparator.comparingInt((node n) -> n.cost)
                            .thenComparing((node n) -> n.state));
        }

        Set<String> explored = new HashSet<>();
        frontier.add(searchTree);

        int nodesExpanded = 0;

        while (!frontier.isEmpty()) {
            node currentNode = frontier.poll();
            nodesExpanded++;

            // Check if the current node is the goal
            if (currentNode.isGoal()) {
                String plan = generatePlan(currentNode);
                int pathCost = calculatePathCost(currentNode);
                return plan + ";" + pathCost + ";" + nodesExpanded;
            }

            explored.add(currentNode.state);

            // Expand the node's children
            for (node child : currentNode.children) {
                if (!explored.contains(child.state)) {
                    frontier.add(child);
                }
            }
        }
        return "NOSOLUTION"; // If no solution found
    }

}
