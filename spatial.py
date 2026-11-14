"""Spatial agreement and tile deduplication used by the GFP detector."""
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

def consensus(predictions,um,minimum=2,radius_um=8.):
    records=sorted((float(y),float(x),int(seed)) for seed,p in predictions.items() for y,x in p)
    if not records:return np.empty((0,2)),[]
    points=np.array([r[:2] for r in records]);seeds=np.array([r[2] for r in records]);distance=cdist(points,points)*um
    groups={i:{i} for i in range(len(points))};owner=np.arange(len(points))
    edges=sorted((distance[i,j],i,j) for i in range(len(points)) for j in range(i+1,len(points)) if seeds[i]!=seeds[j] and distance[i,j]<=radius_um)
    for _,i,j in edges:
        a,b=int(owner[i]),int(owner[j])
        if a==b:continue
        aa,bb=groups[a],groups[b]
        if set(seeds[list(aa)])&set(seeds[list(bb)]):continue
        if distance[np.ix_(sorted(aa),sorted(bb))].max()>radius_um:continue
        groups[a]=aa|bb
        for k in bb:owner[k]=a
        del groups[b]
    kept=[];audit=[]
    for group in groups.values():
        ix=sorted(group)
        if len(ix)<minimum:continue
        #keep an actual predicted center, with coordinate-based ties independent of seed names
        medoid=ix[int(np.argmin(distance[np.ix_(ix,ix)].sum(1)))];kept.append(points[medoid])
        audit.append(dict(center=points[medoid].tolist(),members=[dict(seed=int(seeds[k]),center=points[k].tolist()) for k in ix]))
    return np.array(kept).reshape(-1,2),audit

def deduplicate(points,origins,margins,um):
    if not len(points):
        return np.empty(0,int),np.empty(0,int)
    tree=cKDTree(points)
    owner=np.full(len(points),-1,int);kept=[]
    for i in np.lexsort((np.arange(len(points)),-margins)):
        if owner[i]>=0:
            continue
        kept.append(int(i));owner[i]=i
        candidates=tree.query_ball_point(points[i],8/um)
        candidates=sorted(candidates,key=lambda j:(float(np.linalg.norm(points[j]-points[i])),j))
        represented={tuple(origins[i])}
        for j in candidates:
            tile=tuple(origins[j])
            if owner[j]>=0 or tile in represented or np.linalg.norm(points[j]-points[i])*um>=8:
                continue
            owner[j]=i;represented.add(tile)
    return np.asarray(kept,int),owner

def match_points(first, second, um_per_pixel, radius_um=9.):
    first=np.asarray(first,dtype=float).reshape(-1,2)
    second=np.asarray(second,dtype=float).reshape(-1,2)
    graph=csr_matrix(cdist(first,second)*um_per_pixel<=radius_um)
    return maximum_bipartite_matching(graph,perm_type="column")


def within_circle(points, center, um_per_pixel, radius_um=250.):
    points=np.asarray(points,dtype=float).reshape(-1,2)
    return np.flatnonzero(np.linalg.norm(points-np.asarray(center),axis=1)*um_per_pixel<=radius_um)
