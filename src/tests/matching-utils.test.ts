import { test, describe } from 'node:test';
import assert from 'node:assert/strict';
import { findPerfectMatching } from '../utils/matching-utils.js';

function hasPerfectMatching(graph: number[][]): boolean {
  const mate = new Array<number>(graph.length).fill(-1);
  const search = (): boolean => {
    const i = mate.indexOf(-1);
    if (i < 0) return true;
    for (const j of graph[i]) {
      if (mate[j] < 0 && j !== i) {
        mate[i] = j;
        mate[j] = i;
        if (search()) return true;
        mate[i] = mate[j] = -1;
      }
    }
    return false;
  };
  return search();
}

function isValidPerfectMatching(graph: number[][], matching: number[]): boolean {
  return matching.length === graph.length &&
    matching.every((j, i) => matching[j] === i && graph[i].includes(j));
}

describe('Perfect matching', () => {
  test('graphs with odd cycles (blossoms)', () => {
    // The simple augmenting path search used by the reference implementation
    // returns an invalid matching for this graph
    const graph = [[1, 2, 3], [0, 5, 6, 7], [0, 6, 7], [0, 4, 5], [3, 5], [1, 3, 4], [1, 2], [1, 2]];
    const matching = findPerfectMatching(graph);
    assert(matching !== null);
    assert(isValidPerfectMatching(graph, matching));

    // ...and claims a perfect matching exists for this one, when none does
    const unmatchable = [[3, 4], [2, 6], [1, 4, 5, 7], [0, 4, 6], [0, 2, 3], [2, 6], [1, 3, 5, 7], [2, 6]];
    assert.equal(findPerfectMatching(unmatchable), null);
  });

  test('random graphs agree with brute force', () => {
    let seed = 12345;
    const random = (): number => {
      seed = (seed * 1103515245 + 12345) % 2147483648;
      return seed / 2147483648;
    };

    for (let trial = 0; trial < 5000; trial++) {
      const n = 2 * (1 + Math.floor(random() * 5));
      const p = random() * 0.6;
      const graph: number[][] = Array.from({ length: n }, () => []);
      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          if (random() < p) {
            graph[i].push(j);
            graph[j].push(i);
          }
        }
      }

      const matching = findPerfectMatching(graph);
      assert.equal(matching !== null, hasPerfectMatching(graph), JSON.stringify(graph));
      if (matching !== null) {
        assert(isValidPerfectMatching(graph, matching), JSON.stringify(graph));
      }
    }
  });
});
