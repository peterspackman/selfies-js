/**
 * Molecular graph implementation for SELFIES
 */

import { getBondingCapacity } from './bond-constraints.js';
import { findPerfectMatching } from './utils/matching-utils.js';
import { Attribution } from './types.js';
import type { AttributionMap } from './types.js';
import { AROMATIC_VALENCES, VALENCE_ELECTRONS } from './constants.js';

export { Attribution };

export class Atom {
  /**
   * An atom with associated specifications (e.g. charge, chirality).
   */
  
  public index: number | null = null;
  public element: string;
  public isAromatic: boolean;
  public isotope: number | null;
  public chirality: string | null;
  public hCount: number | null;
  public charge: number;
  
  private _bondingCapacityCache: number | null = null;

  constructor(
    element: string,
    isAromatic: boolean,
    isotope: number | null = null,
    chirality: string | null = null,
    hCount: number | null = null,
    charge: number = 0
  ) {
    this.element = element;
    this.isAromatic = isAromatic;
    this.isotope = isotope;
    this.chirality = chirality;
    this.hCount = hCount;
    this.charge = charge;
  }

  get bondingCapacity(): number {
    if (this._bondingCapacityCache === null) {
      let bondCap = getBondingCapacity(this.element, this.charge);
      bondCap -= (this.hCount === null) ? 0 : this.hCount;
      this._bondingCapacityCache = bondCap;
    }
    return this._bondingCapacityCache;
  }

  clearBondingCapacityCache(): void {
    this._bondingCapacityCache = null;
  }

  invertChirality(): void {
    if (this.chirality === '@') {
      this.chirality = '@@';
    } else if (this.chirality === '@@') {
      this.chirality = '@';
    }
  }
}

export class DirectedBond {
  /**
   * A bond that contains directional information.
   */
  
  public src: number;
  public dst: number;
  public order: number;
  public stereo: string | null;
  public ringBond: boolean;

  constructor(
    src: number,
    dst: number,
    order: number,
    stereo: string | null,
    ringBond: boolean
  ) {
    this.src = src;
    this.dst = dst;
    this.order = order;
    this.stereo = stereo;
    this.ringBond = ringBond;
  }
}

export class MolecularGraph {
  /**
   * A molecular graph.
   * 
   * Molecules can be viewed as weighted undirected graphs. However, SMILES
   * and SELFIES strings are more naturally represented as weighted directed
   * graphs, where the direction of the edges specifies the order of atoms
   * and bonds in the string.
   */
  
  private _roots: number[] = [];
  private _atoms: Atom[] = [];
  private _bondDict: Map<string, DirectedBond> = new Map();
  private _adjList: DirectedBond[][] = [];
  private _bondCounts: number[] = [];
  private _ringBondFlags: boolean[] = [];
  private _delocalSubgraph: Map<number, number[]> = new Map();
  private _attribution: Map<Atom | DirectedBond, Attribution[]> = new Map();
  private _attributable: boolean;

  constructor(attributable: boolean = false) {
    this._attributable = attributable;
  }

  get length(): number {
    return this._atoms.length;
  }

  private _bondKey(a: number, b: number): string {
    return `${a},${b}`;
  }

  hasBond(a: number, b: number): boolean {
    if (a > b) [a, b] = [b, a];
    return this._bondDict.has(this._bondKey(a, b));
  }

  hasOutRingBond(src: number): boolean {
    return this._ringBondFlags[src] || false;
  }

  getAttribution(o: DirectedBond | Atom): Attribution[] | null {
    if (this._attributable && this._attribution.has(o)) {
      return this._attribution.get(o)!;
    }
    return null;
  }

  getRoots(): number[] {
    return [...this._roots];
  }

  getAtom(idx: number): Atom {
    return this._atoms[idx];
  }

  getAtoms(): Atom[] {
    return [...this._atoms];
  }

  getOutDirBonds(src: number): DirectedBond[] {
    return [...this._adjList[src]];
  }

  getDirBond(src: number, dst: number): DirectedBond {
    const key = this._bondKey(src, dst);
    const bond = this._bondDict.get(key);
    if (!bond) {
      throw new Error(`No bond found between ${src} and ${dst}`);
    }
    // Return directed bond in the requested direction
    if (bond.src === src) {
      return bond;
    } else {
      // Create reversed bond
      return new DirectedBond(src, dst, bond.order, bond.stereo, bond.ringBond);
    }
  }

  getBondCount(idx: number): number {
    return this._bondCounts[idx] || 0;
  }

  addAtom(atom: Atom, markRoot: boolean = false): Atom {
    atom.index = this.length;

    if (markRoot) {
      this._roots.push(atom.index);
    }
    
    this._atoms.push(atom);
    this._adjList.push([]);
    this._bondCounts.push(0);
    this._ringBondFlags.push(false);
    
    if (atom.isAromatic) {
      this._delocalSubgraph.set(atom.index, []);
    }
    
    return atom;
  }

  addAttribution(o: DirectedBond | Atom, attr: Attribution[] | null): void {
    if (this._attributable && attr) {
      if (this._attribution.has(o)) {
        this._attribution.get(o)!.push(...attr);
      } else {
        this._attribution.set(o, [...attr]);
      }
    }
  }

  addBond(src: number, dst: number, order: number, stereo: string | null): DirectedBond {
    if (src >= dst) {
      throw new Error(`Source must be less than destination: ${src} >= ${dst}`);
    }

    const bond = new DirectedBond(src, dst, order, stereo, false);
    this._addBondAtLoc(bond, -1);
    this._bondCounts[src] += order;
    this._bondCounts[dst] += order;

    if (order === 1.5) {
      if (!this._delocalSubgraph.has(src)) this._delocalSubgraph.set(src, []);
      if (!this._delocalSubgraph.has(dst)) this._delocalSubgraph.set(dst, []);
      this._delocalSubgraph.get(src)!.push(dst);
      this._delocalSubgraph.get(dst)!.push(src);
    }
    
    return bond;
  }

  addPlaceholderBond(src: number): number {
    const outEdges = this._adjList[src];
    outEdges.push(null as any);
    return outEdges.length - 1;
  }

  addRingBond(
    a: number, 
    b: number,
    order: number,
    aStereo: string | null, 
    bStereo: string | null,
    aPos: number = -1, 
    bPos: number = -1
  ): void {
    const aBond = new DirectedBond(a, b, order, aStereo, true);
    const bBond = new DirectedBond(b, a, order, bStereo, true);
    
    this._addBondAtLoc(aBond, aPos);
    this._addBondAtLoc(bBond, bPos);
    
    this._bondCounts[a] += order;
    this._bondCounts[b] += order;
    this._ringBondFlags[a] = true;
    this._ringBondFlags[b] = true;

    if (order === 1.5) {
      if (!this._delocalSubgraph.has(a)) this._delocalSubgraph.set(a, []);
      if (!this._delocalSubgraph.has(b)) this._delocalSubgraph.set(b, []);
      this._delocalSubgraph.get(a)!.push(b);
      this._delocalSubgraph.get(b)!.push(a);
    }
  }

  updateBondOrder(a: number, b: number, newOrder: number): void {
    if (!(1 <= newOrder && newOrder <= 3)) {
      throw new Error(`Invalid bond order: ${newOrder}`);
    }

    if (a > b) [a, b] = [b, a];
    
    const aToB = this._bondDict.get(this._bondKey(a, b));
    if (!aToB) {
      throw new Error(`Bond not found: ${a} -> ${b}`);
    }
    
    if (newOrder === aToB.order) {
      return;
    }

    let bonds: DirectedBond[];
    if (aToB.ringBond) {
      const bToA = this._bondDict.get(this._bondKey(b, a));
      if (!bToA) {
        throw new Error(`Reverse bond not found: ${b} -> ${a}`);
      }
      bonds = [aToB, bToA];
    } else {
      bonds = [aToB];
    }

    const oldOrder = bonds[0].order;
    for (const bond of bonds) {
      bond.order = newOrder;
    }
    
    this._bondCounts[a] += (newOrder - oldOrder);
    this._bondCounts[b] += (newOrder - oldOrder);
  }

  private _addBondAtLoc(bond: DirectedBond, pos: number): void {
    this._bondDict.set(this._bondKey(bond.src, bond.dst), bond);

    const outEdges = this._adjList[bond.src];
    if (pos === -1 || pos === outEdges.length) {
      outEdges.push(bond);
    } else if (outEdges[pos] === null) {
      outEdges[pos] = bond;
    } else {
      outEdges.splice(pos, 0, bond);
    }
  }

  isKekulized(): boolean {
    return this._delocalSubgraph.size === 0;
  }

  kekulize(): boolean {
    // Algorithm based on Depth-First article by Richard L. Apodaca
    // Reference: https://depth-first.com/articles/2020/02/10/
    // a-comprehensive-treatment-of-aromaticity-in-the-smiles-language/

    if (this.isKekulized()) {
      return true;
    }

    const ds = this._delocalSubgraph;
    const keptNodes = new Set<number>();
    
    // Prune nodes that can't be in a perfect matching
    for (const node of ds.keys()) {
      if (!this._pruneFromDs(node)) {
        keptNodes.add(node);
      }
    }

    // Relabel kept DS nodes to be 0, 1, 2, ...
    const labelToNode = [...keptNodes].sort((a, b) => a - b);
    const nodeToLabel = new Map<number, number>();
    labelToNode.forEach((node, i) => nodeToLabel.set(node, i));

    // Create pruned and relabelled DS
    const prunedDs: number[][] = Array(keptNodes.size).fill(null).map(() => []);
    for (const node of keptNodes) {
      const label = nodeToLabel.get(node)!;
      for (const adj of ds.get(node) || []) {
        if (keptNodes.has(adj)) {
          prunedDs[label].push(nodeToLabel.get(adj)!);
        }
      }
    }

    const matching = findPerfectMatching(prunedDs);
    if (matching === null) {
      return false;
    }

    // De-aromatize and then make double bonds
    for (const node of ds.keys()) {
      for (const adj of ds.get(node) || []) {
        this.updateBondOrder(node, adj, 1);
      }
      if (node < this._atoms.length) {
        this._atoms[node].isAromatic = false;
        this._bondCounts[node] = Math.round(this._bondCounts[node]);
      }
    }

    // Create double bonds based on perfect matching
    for (let i = 0; i < matching.length; i++) {
      const matchedLabel = matching[i];
      if (matchedLabel !== null && i < matchedLabel) { // Process each pair only once
        const node1 = labelToNode[i];
        const node2 = labelToNode[matchedLabel];
        this.updateBondOrder(node1, node2, 2);
      }
    }

    this._delocalSubgraph.clear(); // Clear DS
    return true;
  }

  private _pruneFromDs(node: number): boolean {
    const adjNodes = this._delocalSubgraph.get(node) || [];
    if (adjNodes.length === 0) {
      return true; // Prune isolated nodes
    }
    
    if (node >= this._atoms.length) {
      return true; // Invalid node
    }
    
    const atom = this._atoms[node];
    
    // Check if element can be aromatic
    const aromaticValences = AROMATIC_VALENCES[atom.element];
    if (!aromaticValences) {
      return true; // Element cannot be aromatic
    }
    
    // Calculate electron contribution for aromatic system
    const valenceElectrons = VALENCE_ELECTRONS[atom.element] || 0;
    const charge = atom.charge;
    const hCount = atom.hCount || 0;
    
    // Count non-aromatic bonds (bonds to non-aromatic atoms)
    let usedElectrons = 0;
    let aromaticBonds = 0;
    
    for (let i = 0; i < this._adjList[node].length; i++) {
      const bond = this._adjList[node][i];
      if (bond && bond.dst < this._atoms.length) {
        const dstAtom = this._atoms[bond.dst];
        if (dstAtom.isAromatic && this._delocalSubgraph.has(bond.dst)) {
          aromaticBonds++;
        } else {
          // Non-aromatic bond - count electrons used
          usedElectrons += bond.order;
        }
      }
    }
    
    // Add electrons used for hydrogen atoms
    usedElectrons += hCount;
    
    // Calculate available electrons for aromatic system
    const availableElectrons = valenceElectrons - charge - usedElectrons;
    
    // Each aromatic bond uses 1.5 electrons on average, but we need
    // to check if the atom can contribute the right number of electrons
    // to satisfy one of its allowed aromatic valences
    
    // Check if any aromatic valence is achievable
    let canBeAromatic = false;
    for (const aromaticValence of aromaticValences) {
      const requiredElectrons = aromaticValence - usedElectrons;
      
      // For aromatic systems, we need to check if the electron count works
      // The atom needs to be able to form the aromatic bonds plus contribute
      // to the pi system
      if (requiredElectrons >= aromaticBonds && 
          availableElectrons >= aromaticBonds &&
          availableElectrons <= requiredElectrons) {
        canBeAromatic = true;
        break;
      }
    }
    
    // Special case for common aromatic atoms with typical patterns
    if (!canBeAromatic && aromaticBonds > 0) {
      // Carbon in aromatic rings: typically contributes 1 electron to pi system
      if (atom.element === 'C' && aromaticBonds === 2 && availableElectrons >= 2) {
        canBeAromatic = true;
      }
      // Nitrogen: can contribute 1 or 2 electrons depending on hybridization
      else if (atom.element === 'N' && aromaticBonds <= 2 && availableElectrons >= 1) {
        canBeAromatic = true;
      }
      // Oxygen and Sulfur: typically contribute 2 electrons (lone pair)
      else if ((atom.element === 'O' || atom.element === 'S') && 
               aromaticBonds <= 2 && availableElectrons >= 2) {
        canBeAromatic = true;
      }
    }
    
    return !canBeAromatic; // Prune if cannot be aromatic
  }
}