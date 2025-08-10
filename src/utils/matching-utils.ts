/**
 * Perfect matching algorithm for kekulization
 * Based on the Python SELFIES implementation
 */

export function findPerfectMatching(graph: number[][]): (number | null)[] | null {
  /**
   * Finds a perfect matching for an undirected graph (without self-loops).
   * 
   * @param graph - an adjacency list representing the input graph
   * @returns a list representing a perfect matching, where j is the i-th
   *   element if nodes i and j are matched. Returns null if the graph cannot
   *   be perfectly matched.
   */

  // Start with a maximal matching for efficiency
  const matching = greedyMatching(graph);

  const unmatched = new Set<number>();
  for (let i = 0; i < graph.length; i++) {
    if (matching[i] === null) {
      unmatched.add(i);
    }
  }

  while (unmatched.size > 0) {
    // Find augmenting path which starts at root
    const root = unmatched.values().next().value as number;
    unmatched.delete(root);
    
    const path = findAugmentingPath(graph, root, matching);

    if (path === null) {
      return null;
    } else {
      flipAugmentingPath(matching, path);
      unmatched.delete(path[0]);
      unmatched.delete(path[path.length - 1]);
    }
  }

  return matching;
}

function greedyMatching(graph: number[][]): (number | null)[] {
  const matching: (number | null)[] = new Array(graph.length).fill(null);
  const freeDegrees = graph.map(neighbors => neighbors.length);

  // Create priority queue of nodes (prioritize nodes with fewer unmatched neighbors)
  const nodePQueue: Array<{degree: number, node: number}> = [];
  for (let i = 0; i < graph.length; i++) {
    nodePQueue.push({degree: freeDegrees[i], node: i});
  }
  nodePQueue.sort((a, b) => a.degree - b.degree);

  let queueIndex = 0;
  while (queueIndex < nodePQueue.length) {
    const {node} = nodePQueue[queueIndex++];

    if (matching[node] !== null || freeDegrees[node] === 0) {
      continue; // node cannot be matched
    }

    // Find neighbor with smallest free degree
    let bestNeighbor = -1;
    let bestDegree = Infinity;
    
    for (const neighbor of graph[node]) {
      if (matching[neighbor] === null && freeDegrees[neighbor] < bestDegree) {
        bestNeighbor = neighbor;
        bestDegree = freeDegrees[neighbor];
      }
    }

    if (bestNeighbor !== -1) {
      // Match node with best neighbor
      matching[node] = bestNeighbor;
      matching[bestNeighbor] = node;

      // Update free degrees for all neighbors
      for (const neighbor of graph[node]) {
        freeDegrees[neighbor]--;
      }
      for (const neighbor of graph[bestNeighbor]) {
        freeDegrees[neighbor]--;
      }
    }
  }

  return matching;
}

function findAugmentingPath(
  graph: number[][], 
  root: number, 
  matching: (number | null)[]
): number[] | null {
  const visited = new Set<number>();
  const parent: (number | null)[] = new Array(graph.length).fill(null);
  const queue: number[] = [root];
  visited.add(root);

  while (queue.length > 0) {
    const node = queue.shift()!;

    for (const neighbor of graph[node]) {
      if (visited.has(neighbor)) {
        continue;
      }

      visited.add(neighbor);
      parent[neighbor] = node;

      const match = matching[neighbor];
      if (match === null) {
        // Found augmenting path - reconstruct it
        const path: number[] = [neighbor];
        let current = neighbor;
        while (parent[current] !== null) {
          current = parent[current]!;
          path.unshift(current);
        }
        return path;
      } else {
        // Continue BFS through the matched edge
        if (!visited.has(match)) {
          visited.add(match);
          parent[match] = neighbor;
          queue.push(match);
        }
      }
    }
  }

  return null; // No augmenting path found
}

function flipAugmentingPath(matching: (number | null)[], path: number[]): void {
  for (let i = 0; i < path.length - 1; i++) {
    const u = path[i];
    const v = path[i + 1];
    matching[u] = v;
    matching[v] = u;
  }
}