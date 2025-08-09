/**
 * Molecular graph implementation for SELFIES
 */

import { getBondingCapacity } from './bond-constraints.js';
import { Attribution } from './types.js';
import type { AttributionMap } from './types.js';

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

  getDirBond(src: number, dst: number): DirectedBond {
    const bond = this._bondDict.get(this._bondKey(src, dst));
    if (!bond) {
      throw new Error(`Bond not found: ${src} -> ${dst}`);
    }
    return bond;
  }

  getOutDirBonds(src: number): DirectedBond[] {
    return [...this._adjList[src]];
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
    if (this.isKekulized()) {
      return true;
    }

    const visited = new Set<number>();
    
    for (const [atom, neighbors] of this._delocalSubgraph) {
      if (visited.has(atom)) continue;
      
      let bondOrder = 1;
      for (const neighbor of neighbors) {
        if (visited.has(neighbor)) continue;
        
        const key = this._bondKey(Math.min(atom, neighbor), Math.max(atom, neighbor));
        const bond = this._bondDict.get(key);
        if (bond && bond.order === 1.5) {
          bond.order = bondOrder;
          bondOrder = bondOrder === 1 ? 2 : 1;
        }
        
        visited.add(neighbor);
      }
      visited.add(atom);
    }

    this._delocalSubgraph.clear();
    return true;
  }
}