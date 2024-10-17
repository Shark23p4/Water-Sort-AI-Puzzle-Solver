
import java.util.*;

class node {

    String state; // Representing the state of the node
    List<node> children; // List of child nodes
    boolean isGoal; // Indicator for goal state
    String action; // Action that led to this node (for plan generation)
    node parent; // Reference to the parent node for path reconstruction

    // Constructor
    public node(String state) {
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

class WaterSortSearch {

    List<Stack<Character>> bottles;
    int bottleCapacity = 0;

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
        node root = new node(encodeState(bottles));
        Set<String> visitedStates = new HashSet<>();
        generateMoves(root, bottles, visitedStates);
        return root;
    }

    // Recursive method to generate possible moves for a given state, avoiding
    // repeated states and stopping at goal state
    private void generateMoves(node node, List<Stack<Character>> currBottles,
            Set<String> visitedStates) {
        int numberOfBottles = bottles.size();

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
                if (i != j && validMove(bottles, i, j, bottleCapacity)) {
                    List<Stack<Character>> newBottles = Pour(i, j);
                    String newState = encodeState(newBottles);

                    // Only proceed if this new state hasn't been visited yet
                    if (!visitedStates.contains(newState)) {
                        node child = new node(newState);
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

    public String solve(String initialState, String strategy, boolean visualize) {
        // Parse user input
        int heuristicType = 0;
        if (strategy.matches("GR[12]") || strategy.matches("AS[12]")) {
            heuristicType = Character.getNumericValue(strategy.charAt(2));
            strategy = strategy.substring(0, 2); // Trim heuristic type
        }

        long startTime = System.nanoTime(); // Start tracking time
        long startMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory(); // Start tracking
        // memory

        node searchTree = generateMoveTree(initialState);
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
        String initialState = "5;4;r,g,e,e;b,g,e,e;y,e,e,e;b,r,e,e;r,g,e,e";

        String searchType = "BF"; // Example input startegy

        boolean visualize = true;

        WaterSortSearch WaterSortSearch = new WaterSortSearch();

        WaterSortSearch.solve(initialState, searchType, visualize);

    }
}

class BFS {

    public static String search(node searchTree, boolean visualize) {
        // BFS Implementation
        int nodesExpanded = 0;
        StringBuilder plan = new StringBuilder();
        int pathCost = 0;

        Queue<TreeNode> frontier = new LinkedList<>();
        Set<TreeNode> explored = new HashSet<>();
        frontier.add(searchTree);

        while (!frontier.isEmpty()) {
            TreeNode currentNode = frontier.poll();
            explored.add(currentNode);
            nodesExpanded++;

            // Visualize current state
            if (visualize) {
                System.out.println("Exploring: " + currentNode.state);
            }

            // Check if the current node is the goal
            if (currentNode.isGoal()) {
                plan.append(generatePlan(currentNode)); // Generate plan from goal node
                pathCost = calculatePathCost(currentNode); // Calculate path cost
                return plan.toString() + ";" + pathCost + ";" + nodesExpanded; // Return the result
            }

            // Add all unvisited children to the frontier
            for (TreeNode neighbor : currentNode.children) {
                if (!explored.contains(neighbor) && !frontier.contains(neighbor)) {
                    neighbor.action = currentNode.action; // Set the action that led to this node
                    frontier.add(neighbor);
                }
            }
        }

        return "NOSOLUTION"; // Return if no solution is found
    }

    private static String generatePlan(TreeNode node) {
        // Generate a solution plan based on the node's path from the root
        List<String> moves = new ArrayList<>();
        TreeNode current = node;

        // Assuming each node has a reference to its parent (not shown in TreeNode)
        while (current != null) {
            if (!current.action.isEmpty()) {
                moves.add(current.action); // Collect the action
            }
            current = null; // Replace with current.parent if parent reference is added
        }

        Collections.reverse(moves); // Reverse to get path from root to goal
        return String.join(",", moves); // Join moves with commas
    }

    private static int calculatePathCost(TreeNode node) {
        // Calculate path cost based on the depth of the node
        int cost = 0;
        TreeNode current = node;

        // Assuming each move has a uniform cost of 1 (or customize if needed)
        while (current != null) {
            cost += 1; // Increment cost for each level
            current = null; // Replace with current.parent if parent reference is added
        }
        return cost - 1; // Subtract 1 to avoid counting the goal node itself
    }

}
