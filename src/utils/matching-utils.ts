/**
 * Perfect matching algorithm for kekulization
 * Based on the Python SELFIES implementation
 */

/**
 * Finds a perfect matching for an undirected graph (without self-loops).
 *
 * @param graph - an adjacency list representing the input graph
 * @returns a list representing a perfect matching, where j is the i-th
 *   element if nodes i and j are matched. Returns null if the graph cannot
 *   be perfectly matched.
 */
export function findPerfectMatching(graph: number[][]): number[] | null {
  // Start with a maximal matching for efficiency
  const greedy = greedyMatching(graph);
  const matching = [...greedy];

  // Visit unmatched nodes in ascending order (as the reference implementation does)
  const unmatched = new Set<number>();
  for (let i = 0; i < graph.length; i++) {
    if (matching[i] === null) {
      unmatched.add(i);
    }
  }

  let ok = true;
  for (const root of unmatched) {
    unmatched.delete(root);
    const path = findAugmentingPath(graph, root, matching);
    if (path === null) {
      ok = false;
      break;
    }
    flipAugmentingPath(matching, path);
    unmatched.delete(path[0]);
    unmatched.delete(path[path.length - 1]);
  }

  if (ok && isPerfectMatching(graph, matching)) {
    return matching as number[];
  }

  // The simple search above does not handle odd cycles (blossoms): it can miss
  // augmenting paths, or flip a path that revisits a node and so corrupt the
  // matching. Fall back to Edmonds' algorithm from the greedy matching.
  return blossomPerfectMatching(graph, greedy);
}

function isPerfectMatching(graph: number[][], matching: (number | null)[]): boolean {
  return matching.every((mate, i) =>
    mate !== null && matching[mate] === i && graph[i].includes(mate));
}

/**
 * Minimal binary heap over [priority, node] pairs, ordered like Python tuples.
 */
class NodeHeap {
  private items: [number, number][] = [];

  get size(): number {
    return this.items.length;
  }

  private less(a: [number, number], b: [number, number]): boolean {
    return a[0] < b[0] || (a[0] === b[0] && a[1] < b[1]);
  }

  push(item: [number, number]): void {
    const items = this.items;
    items.push(item);
    let i = items.length - 1;
    while (i > 0) {
      const parent = (i - 1) >> 1;
      if (!this.less(items[i], items[parent])) break;
      [items[i], items[parent]] = [items[parent], items[i]];
      i = parent;
    }
  }

  pop(): [number, number] {
    const items = this.items;
    const top = items[0];
    const last = items.pop()!;
    if (items.length > 0) {
      items[0] = last;
      let i = 0;
      for (;;) {
        const l = 2 * i + 1;
        const r = l + 1;
        let smallest = i;
        if (l < items.length && this.less(items[l], items[smallest])) smallest = l;
        if (r < items.length && this.less(items[r], items[smallest])) smallest = r;
        if (smallest === i) break;
        [items[i], items[smallest]] = [items[smallest], items[i]];
        i = smallest;
      }
    }
    return top;
  }
}

function greedyMatching(graph: number[][]): (number | null)[] {
  const matching: (number | null)[] = new Array(graph.length).fill(null);
  // freeDegrees[i] = number of unmatched neighbors for node i
  const freeDegrees = graph.map(neighbors => neighbors.length);

  // Prioritize nodes with fewer unmatched neighbors
  const nodePQueue = new NodeHeap();
  for (let i = 0; i < graph.length; i++) {
    nodePQueue.push([freeDegrees[i], i]);
  }

  while (nodePQueue.size > 0) {
    const [, node] = nodePQueue.pop();
    if (matching[node] !== null || freeDegrees[node] === 0) {
      continue; // node cannot be matched
    }

    // Match node with first unmatched neighbor
    const mate = graph[node].find(i => matching[i] === null)!;
    matching[node] = mate;
    matching[mate] = node;

    for (const adj of [...graph[node], ...graph[mate]]) {
      freeDegrees[adj] -= 1;
      if (matching[adj] === null && freeDegrees[adj] > 0) {
        nodePQueue.push([freeDegrees[adj], adj]);
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
  // Run modified BFS to find path from root to unmatched node
  let otherEnd: number | null = null;
  const nodeQueue: number[] = [root];
  let head = 0;

  // Parent BFS tree - null indicates an unvisited node
  const parents: ([number, number] | null)[] = new Array(graph.length).fill(null);
  parents[root] = [-1, -1];

  while (head < nodeQueue.length) {
    const node = nodeQueue[head++];
    for (const adj of graph[node]) {
      const adjMate = matching[adj];
      if (adjMate === null) {
        if (adj !== root) { // augmenting path found!
          parents[adj] = [node, adj];
          otherEnd = adj;
          break;
        }
      } else if (parents[adjMate] === null) { // adjMate not visited
        parents[adjMate] = [node, adj];
        nodeQueue.push(adjMate);
      }
    }
    if (otherEnd !== null) {
      break;
    }
  }

  if (otherEnd === null) {
    return null;
  }

  const path: number[] = [];
  let node = otherEnd;
  while (node !== root) {
    const [prev, via] = parents[node]!;
    path.push(via, prev);
    node = prev;
  }
  return path;
}

function flipAugmentingPath(matching: (number | null)[], path: number[]): void {
  for (let i = 0; i < path.length; i += 2) {
    const a = path[i];
    const b = path[i + 1];
    matching[a] = b;
    matching[b] = a;
  }
}

/**
 * Extends a partial matching to a perfect matching using Edmonds' blossom
 * algorithm, or returns null if no perfect matching exists.
 */
function blossomPerfectMatching(
  graph: number[][],
  initial: (number | null)[]
): number[] | null {
  const n = graph.length;
  const match = initial.map(m => (m === null ? -1 : m));
  const parent = new Array<number>(n);
  const base = new Array<number>(n);
  const used = new Array<boolean>(n);
  const blossom = new Array<boolean>(n);

  const lca = (a: number, b: number): number => {
    const seen = new Array<boolean>(n).fill(false);
    for (;;) {
      a = base[a];
      seen[a] = true;
      if (match[a] === -1) break;
      a = parent[match[a]];
    }
    for (;;) {
      b = base[b];
      if (seen[b]) return b;
      b = parent[match[b]];
    }
  };

  const markPath = (v: number, b: number, child: number): void => {
    while (base[v] !== b) {
      blossom[base[v]] = true;
      blossom[base[match[v]]] = true;
      parent[v] = child;
      child = match[v];
      v = parent[match[v]];
    }
  };

  const findPath = (root: number): number => {
    used.fill(false);
    parent.fill(-1);
    for (let i = 0; i < n; i++) base[i] = i;
    used[root] = true;
    const queue = [root];
    let head = 0;
    while (head < queue.length) {
      const v = queue[head++];
      for (const to of graph[v]) {
        if (base[v] === base[to] || match[v] === to) continue;
        if (to === root || (match[to] !== -1 && parent[match[to]] !== -1)) {
          // Found a blossom; contract it
          const curBase = lca(v, to);
          blossom.fill(false);
          markPath(v, curBase, to);
          markPath(to, curBase, v);
          for (let i = 0; i < n; i++) {
            if (blossom[base[i]]) {
              base[i] = curBase;
              if (!used[i]) {
                used[i] = true;
                queue.push(i);
              }
            }
          }
        } else if (parent[to] === -1) {
          parent[to] = v;
          if (match[to] === -1) return to;
          used[match[to]] = true;
          queue.push(match[to]);
        }
      }
    }
    return -1;
  };

  for (let root = 0; root < n; root++) {
    if (match[root] !== -1) continue;
    let v = findPath(root);
    if (v === -1) return null;
    while (v !== -1) {
      const pv = parent[v];
      const ppv = match[pv];
      match[v] = pv;
      match[pv] = v;
      v = ppv;
    }
  }

  return match;
}
