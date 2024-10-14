import java.util.*;

class TreeNode {
    String state; // Representing the state of the node
    List<TreeNode> children; // List of child nodes
    boolean isGoal; // Indicator for goal state
    String action; // Action that led to this node (for plan generation)
    TreeNode parent; // Reference to the parent node for path reconstruction

    // Constructor
    public TreeNode(String state) {
        this.state = state;
        this.children = new ArrayList<>();
        this.isGoal = false; // Default to false
        this.action = ""; // Initialize with no action
        this.parent = null; // Initialize with no parent
    }

    // Method to set the goal status
    public void setGoal(boolean goalStatus) {
        this.isGoal = goalStatus;
    }

    public boolean isGoal() {
        return this.isGoal;
    }
}

class WaterSortPuzzleSolver {
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
    public static boolean validMove(List<Stack<Character>> bottles, int i, int j, int bottleCapacity) {
        if (bottles.get(i).isEmpty())
            return false; // Bottle i is empty
        if (bottles.get(j).size() == bottleCapacity)
            return false; // Bottle j is full
        if (bottles.get(j).isEmpty() || bottles.get(i).peek().equals(bottles.get(j).peek()))
            return true;
        return false;
    }

    // Method to perform a move and return the new state as a string
    public static List<Stack<Character>> performMove(List<Stack<Character>> bottles, int i, int j, int bottleCapacity) {
        List<Stack<Character>> newBottles = deepCopy(bottles);
        char color = newBottles.get(i).peek();
        int layersToMove = 0;

        // Count how many contiguous layers of the same color are on top of bottle i
        while (layersToMove < newBottles.get(i).size()
                && newBottles.get(i).get(newBottles.get(i).size() - 1 - layersToMove) == color) {
            layersToMove++;
        }

        int emptySpaceInJ = bottleCapacity - newBottles.get(j).size();
        int moveLayers = Math.min(layersToMove, emptySpaceInJ);

        // Perform the move by transferring the layers
        for (int k = 0; k < moveLayers; k++) {
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
    public static String encodeState(List<Stack<Character>> bottles, int bottleCapacity) {
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
    public static TreeNode generateMoveTree(String initialState) {
        List<Stack<Character>> initialBottles = parseInitialState(initialState);
        int bottleCapacity = getBottleCapacity(initialState); // Derive bottleCapacity from initialState
        TreeNode root = new TreeNode(encodeState(initialBottles, bottleCapacity));
        Set<String> visitedStates = new HashSet<>();
        generateMoves(root, initialBottles, bottleCapacity, visitedStates);
        return root;
    }

    // Recursive method to generate possible moves for a given state, avoiding
    // repeated states and stopping at goal state
    private static void generateMoves(TreeNode node, List<Stack<Character>> bottles, int bottleCapacity,
            Set<String> visitedStates) {
        int numberOfBottles = bottles.size();

        // Check if the current state is a goal state
        if (isGoalState(bottles)) {
            node.setGoal(true);
            return; // Stop recursion if the goal state is reached
        }

        String currentState = encodeState(bottles, bottleCapacity);
        visitedStates.add(currentState); // Mark current state as visited

        // Generate all valid moves and create child nodes
        for (int i = 0; i < numberOfBottles; i++) {
            for (int j = 0; j < numberOfBottles; j++) {
                if (i != j && validMove(bottles, i, j, bottleCapacity)) {
                    List<Stack<Character>> newBottles = performMove(bottles, i, j, bottleCapacity);
                    String newState = encodeState(newBottles, bottleCapacity);

                    // Only proceed if this new state hasn't been visited yet
                    if (!visitedStates.contains(newState)) {
                        TreeNode child = new TreeNode(newState);
                        child.parent = node;
                        child.action = "Pour_" + i + "_" + j;
                        node.children.add(child);

                        // Recurse for further moves
                        generateMoves(child, newBottles, bottleCapacity, visitedStates);
                    }
                }
            }
        }
    }

    public static String solve(String initialState, String strategy, boolean visualize) {
        // Parse user input
        int heuristicType = 0;
        if (strategy.matches("GR[12]") || strategy.matches("AS[12]")) {
            heuristicType = Character.getNumericValue(strategy.charAt(2));
            strategy = strategy.substring(0, 2); // Trim heuristic type
        }

        long startTime = System.nanoTime(); // Start tracking time
        long startMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory(); // Start tracking
                                                                                                   // memory

        TreeNode searchTree = generateMoveTree(initialState);
        String solution = "";
        switch (strategy) {
            case "BF":
                solution = BFS.search(searchTree, visualize);
                break;
            case "DF":
                // solution = DFS.search(searchTree, visualize);
                break;
            case "ID":
                // solution = IDS.search(searchTree, visualize);
                break;
            case "UC":
                // solution = UCS.search(searchTree, visualize);
                break;
            case "GR":
                // solution = GreedySearch.search(searchTree, heuristicType, visualize);
                break;
            case "AS":
                // solution = AStarSearch.search(searchTree, heuristicType, visualize);
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
            if(solutionParts.length > 1) {
            System.out.println("pathCost: " + solutionParts[1]);
            System.out.println("nodesExpanded: " + solutionParts[2]);
        }
    }

        return solution;
    }

    public static void printTree(TreeNode node, int depth) {
        if (node == null) {
            return;
        }

        // Print the current node's state and moves
        System.out.println("Depth " + depth + ": " + node.state);
        for (TreeNode child : node.children) {
            System.out.println("    Move: " + child.move + ", Resulting State: " + child.state);
        }

        // Recursively print children
        for (TreeNode child : node.children) {
            printTree(child, depth + 1);
        }
    }

    public static void main(String[] args) {
        String initialState = "5;4;r,g,g,g;b,g,g,g;y,y,r,g;g,g,g,g;r,g,r,e,e";

        String searchType = "BF"; // Example input startegy

        boolean visualize = true;

        solve(initialState, searchType, visualize);

    }
}

class BFS {
    public static String search(TreeNode searchTree, boolean visualize) {
        // BFS Implementation
        int nodesExpanded = 0;
        StringBuilder plan = new StringBuilder();
        int pathCost = 0;

        Queue<TreeNode> frontier = new LinkedList<>();
        frontier.add(searchTree);

        while (!frontier.isEmpty()) {
            TreeNode currentNode = frontier.poll();
            nodesExpanded++;

            // Check if the current node is the goal
            if (currentNode.isGoal()) {
                plan.append(generatePlan(currentNode)); // Generate plan from goal node
                pathCost = calculatePathCost(currentNode); // Calculate path cost
                return plan.toString() + ";" + pathCost + ";" + nodesExpanded; // Return the result
            }

            // Add all children to the frontier
            for (TreeNode neighbor : currentNode.children) {
                neighbor.parent = currentNode; // Set the parent for path reconstruction
                frontier.add(neighbor);
            }
        }

        return "NOSOLUTION"; // Return if no solution is found
    }

    private static String generatePlan(TreeNode node) {
        // Generate a solution plan based on the node's path from the root
        List<String> moves = new ArrayList<>();
        TreeNode current = node;

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

    private static int calculatePathCost(TreeNode node) {
        // Calculate path cost based on the depth of the node
        int cost = 0;
        TreeNode current = node;

        // Count the number of edges from the root to the goal node
        while (current != null) {
            cost += 1; // Increment cost for each level
            current = current.parent; // Move to the parent node
        }
        return cost - 1; // Subtract 1 to avoid counting the goal node itself
    }
}

/*
 * // Implementing DFS
 * class DFS {
 * public static String search(String initialState, boolean visualize) {
 * // DFS Implementation
 * int nodesExpanded = 0;
 * String plan = "";
 * int pathCost = 0;
 * 
 * Stack<String> frontier = new Stack<>();
 * Set<String> explored = new HashSet<>();
 * frontier.push(initialState);
 * 
 * while (!frontier.isEmpty()) {
 * String state = frontier.pop();
 * explored.add(state);
 * nodesExpanded++;
 * 
 * if (isGoal(state)) {
 * plan = generatePlan(state);
 * pathCost = calculatePathCost(state);
 * break;
 * }
 * 
 * for (String neighbor : getNeighbors(state)) {
 * if (!explored.contains(neighbor) && !frontier.contains(neighbor)) {
 * frontier.push(neighbor);
 * }
 * }
 * }
 * 
 * if (visualize) {
 * System.out.println("DFS Plan: " + plan);
 * }
 * 
 * return plan + ";" + pathCost + ";" + nodesExpanded;
 * }
 * 
 * private static boolean isGoal(String state) {
 * return state.equals("goal");
 * }
 * 
 * private static String generatePlan(String state) {
 * return "move1, move2, move3";
 * }
 * 
 * private static int calculatePathCost(String state) {
 * return 15;
 * }
 * 
 * private static List<String> getNeighbors(String state) {
 * return Arrays.asList("neighbor1", "neighbor2");
 * }
 * }
 * 
 * // Implementing IDS (Iterative Deepening Search)
 * class IDS {
 * public static String search(String initialState, boolean visualize) {
 * // IDS Implementation
 * int depth = 0;
 * int nodesExpanded = 0;
 * String plan = "";
 * int pathCost = 0;
 * 
 * while (true) {
 * String result = depthLimitedSearch(initialState, depth);
 * if (result != null) {
 * plan = result;
 * pathCost = calculatePathCost(result);
 * break;
 * }
 * depth++;
 * }
 * 
 * nodesExpanded += depth; // Dummy node expansion count
 * 
 * if (visualize) {
 * System.out.println("IDS Plan: " + plan);
 * }
 * 
 * return plan + ";" + pathCost + ";" + nodesExpanded;
 * }
 * 
 * private static String depthLimitedSearch(String state, int depth) {
 * // Perform depth-limited search
 * return state.equals("goal") ? "move1, move2, move3" : null;
 * }
 * 
 * private static int calculatePathCost(String state) {
 * return 12;
 * }
 * }
 * 
 * // Implementing UCS (Uniform Cost Search)
 * class UCS {
 * public static String search(String initialState, boolean visualize) {
 * // UCS Implementation
 * int nodesExpanded = 0;
 * String plan = "";
 * int pathCost = 0;
 * 
 * PriorityQueue<Node> frontier = new PriorityQueue<>(Comparator.comparingInt(n
 * -> n.pathCost));
 * Set<String> explored = new HashSet<>();
 * frontier.add(new Node(initialState, 0));
 * 
 * while (!frontier.isEmpty()) {
 * Node currentNode = frontier.poll();
 * String state = currentNode.state;
 * pathCost = currentNode.pathCost;
 * explored.add(state);
 * nodesExpanded++;
 * 
 * if (isGoal(state)) {
 * plan = generatePlan(state);
 * break;
 * }
 * 
 * for (String neighbor : getNeighbors(state)) {
 * if (!explored.contains(neighbor)) {
 * frontier.add(new Node(neighbor, pathCost + 1)); // Assuming each move costs 1
 * }
 * }
 * }
 * 
 * if (visualize) {
 * System.out.println("UCS Plan: " + plan);
 * }
 * 
 * return plan + ";" + pathCost + ";" + nodesExpanded;
 * }
 * 
 * private static boolean isGoal(String state) {
 * return state.equals("goal");
 * }
 * 
 * private static String generatePlan(String state) {
 * return "move1, move2, move3";
 * }
 * 
 * private static List<String> getNeighbors(String state) {
 * return Arrays.asList("neighbor1", "neighbor2");
 * }
 * 
 * static class Node {
 * String state;
 * int pathCost;
 * 
 * Node(String state, int pathCost) {
 * this.state = state;
 * this.pathCost = pathCost;
 * }
 * }
 * }
 * 
 * // Implementing Greedy Search
 * class GreedySearch {
 * public static String search(String initialState, int heuristicType, boolean
 * visualize) {
 * // Greedy Search Implementation
 * int nodesExpanded = 0;
 * String plan = "";
 * int pathCost = 0;
 * 
 * PriorityQueue<Node> frontier = new PriorityQueue<>(
 * Comparator.comparingInt(n -> heuristic(n.state, heuristicType)));
 * Set<String> explored = new HashSet<>();
 * frontier.add(new Node(initialState, 0));
 * 
 * while (!frontier.isEmpty()) {
 * Node currentNode = frontier.poll();
 * String state = currentNode.state;
 * pathCost = currentNode.pathCost;
 * explored.add(state);
 * nodesExpanded++;
 * 
 * if (isGoal(state)) {
 * plan = generatePlan(state);
 * break;
 * }
 * 
 * for (String neighbor : getNeighbors(state)) {
 * if (!explored.contains(neighbor)) {
 * frontier.add(new Node(neighbor, pathCost + 1));
 * }
 * }
 * }
 * 
 * if (visualize) {
 * System.out.println("Greedy Search Plan: " + plan);
 * }
 * 
 * return plan + ";" + pathCost + ";" + nodesExpanded;
 * }
 * 
 * private static int heuristic(String state, int heuristicType) {
 * // Heuristic calculation based on type (dummy example)
 * if (heuristicType == 1) {
 * return state.length(); // Heuristic 1 (dummy)
 * } else {
 * return state.chars().sum(); // Heuristic 2 (dummy)
 * }
 * }
 * 
 * private static boolean isGoal(String state) {
 * return state.equals("goal");
 * }
 * 
 * private static String generatePlan(String state) {
 * return "move1, move2, move3";
 * }
 * 
 * private static List<String> getNeighbors(String state) {
 * return Arrays.asList("neighbor1", "neighbor2");
 * }
 * 
 * static class Node {
 * String state;
 * int pathCost;
 * 
 * Node(String state, int pathCost) {
 * this.state = state;
 * this.pathCost = pathCost;
 * }
 * }
 * }
 * 
 * // Implementing A* Search
 * class AStarSearch {
 * public static String search(String initialState, int heuristicType, boolean
 * visualize) {
 * // A* Search Implementation
 * int nodesExpanded = 0;
 * String plan = "";
 * int pathCost = 0;
 * 
 * PriorityQueue<Node> frontier = new PriorityQueue<>(
 * Comparator.comparingInt(n -> n.pathCost + heuristic(n.state,
 * heuristicType)));
 * Set<String> explored = new HashSet<>();
 * frontier.add(new Node(initialState, 0));
 * 
 * while (!frontier.isEmpty()) {
 * Node currentNode = frontier.poll();
 * String state = currentNode.state;
 * pathCost = currentNode.pathCost;
 * explored.add(state);
 * nodesExpanded++;
 * 
 * if (isGoal(state)) {
 * plan = generatePlan(state);
 * break;
 * }
 * 
 * for (String neighbor : getNeighbors(state)) {
 * if (!explored.contains(neighbor)) {
 * frontier.add(new Node(neighbor, pathCost + 1));
 * }
 * }
 * }
 * 
 * if (visualize) {
 * System.out.println("A* Plan: " + plan);
 * }
 * 
 * return plan + ";" + pathCost + ";" + nodesExpanded;
 * }
 * 
 * private static int heuristic(String state, int heuristicType) {
 * // Heuristic calculation based on type (dummy example)
 * if (heuristicType == 1) {
 * return state.length(); // Heuristic 1 (dummy)
 * } else {
 * return state.chars().sum(); // Heuristic 2 (dummy)
 * }
 * }
 * 
 * private static boolean isGoal(String state) {
 * return state.equals("goal");
 * }
 * 
 * private static String generatePlan(String state) {
 * return "move1, move2, move3";
 * }
 * 
 * private static List<String> getNeighbors(String state) {
 * return Arrays.asList("neighbor1", "neighbor2");
 * }
 * 
 * static class Node {
 * String state;
 * int pathCost;
 * 
 * Node(String state, int pathCost) {
 * this.state = state;
 * this.pathCost = pathCost;
 * }
 * }
 * }
 */