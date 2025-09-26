//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 16:34:29 taubin>
//------------------------------------------------------------------------
//
// PolygonMesh.cpp
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//     * Redistributions of source code must retain the above
//       copyright notice, this list of conditions and the following
//       disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials
//       provided with the distribution.
//     * Neither the name of the Brown University nor the names of its
//       contributors may be used to endorse or promote products
//       derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GABRIEL
// TAUBIN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include <iostream>
#include "PolygonMesh.hpp"
#include "Partition.hpp"

PolygonMesh::PolygonMesh(const int nVertices, const vector<int>& coordIndex):
  HalfEdges(nVertices,coordIndex),
  _nPartsVertex(),
  _isBoundaryVertex()
{
  int nV = getNumberOfVertices();
  int nE = getNumberOfEdges(); // Edges method
  int nC = getNumberOfCorners();
  int iE;

  // 1) classify the vertices as boundary or internal
  int iV;
  for(iV=0;iV<nV;iV++)
    _isBoundaryVertex.push_back(false);

  // - for edge boundary iE label its two end vertices as boundary

  for(iE = 0; iE < nE; ++iE){
        if(getNumberOfEdgeHalfEdges(iE) == 1){
          _isBoundaryVertex[getVertex0(iE)] = true;
          _isBoundaryVertex[getVertex1(iE)] = true;
        }
  }
  
  // 2) create a partition of the corners in the stack
  Partition partition(nC);

  // 3) for each regular edge
  int e1; int e2;
  for(iE = 0; iE < nE; ++iE){
      if(getNumberOfEdgeHalfEdges(iE) == 2){
          //    - get the two half edges incident to the edge
          //    - join the two pairs of corresponding corners accross the edge
          //    - you need to take into account the relative orientation of
          //      the two incident half edges
          e1 = getEdgeHalfEdge(iE, 0);
          e2 = getEdgeHalfEdge(iE, 1);

          // for the moment let's assume that the mesh does not have
          // singular edges, and that pairs of corners corresponding to the
          // same vertex across inconsistently oriented faces will be joined
          partition.join(e1,getNext(e2));
          partition.join(e2,getNext(e1));

          // note that the partition will end up with the corner separators as
          // singletons, but it doesn't matter for the last step, and
          // the partition will be deleteted upon return

      }
  }
  
  // 4) count number of parts per vertex
  //   - initialize _nPartsVertex array to 0's
  _nPartsVertex.resize(nV,0);
  //    - for each corner iC which is a representative of its subset
  int pId;
  vector<bool> visited(partition.getNumberOfParts(), false);
  cout << visited.size() << endl;
  for(int iC = 0; iC < nC; ++iC)
      if(_coordIndex[iC] >= 0 && (pId = partition.find(iC)) != -1)
          if(!visited[pId]){
             // - get the corresponding vertex index iV and increment _nPartsVertex[iV]
             // - note that all the corners in each subset share a common
            //  vertex index, but multiple subsets may correspond to the
            //  same vertex index, indicating that the vertex is singular
             cout << coordIndex[iC] << " " << pId << endl;
            _nPartsVertex[coordIndex[iC]]++;
            visited[pId] = true;
          }
}


int PolygonMesh::getNumberOfFaces() const {
    return _face[_face.size() - 2]  + 1;
}

int PolygonMesh::getNumberOfEdgeFaces(const int iE) const {
  return getNumberOfEdgeHalfEdges(iE);
}


int PolygonMesh::getEdgeFace(const int iE, const int j) const {
    return getFace(getEdgeHalfEdge(iE, j));
}
// if the arguments fall within their respective ranges, this method
// returns returns true if iF is found in the list of faces incident
// to the edge iE; otherwise it returns false


bool PolygonMesh::isEdgeFace(const int iE, const int iF) const {
    if(iF >= 0 && iF < getNumberOfFaces())
        for(u_int i = 0; i < getNumberOfEdgeFaces(iE); ++i)
            if(getEdgeFace(iE, i) == iF)
                return true;
    return false;
}

// classification of edges

bool PolygonMesh::isBoundaryEdge(const int iE) const {
    return isValidEdge(iE) && getNumberOfEdgeHalfEdges(iE) == 1;
}

bool PolygonMesh::isRegularEdge(const int iE) const {
  return isValidEdge(iE) && getNumberOfEdgeHalfEdges(iE) == 2;
}

bool PolygonMesh::isSingularEdge(const int iE) const {
    return isValidEdge(iE) && getNumberOfEdgeHalfEdges(iE) > 2;
}

// classification of vertices

bool PolygonMesh::isBoundaryVertex(const int iV) const {
  int nV = getNumberOfVertices();
  return (0<=iV && iV<nV)?_isBoundaryVertex[iV]:false;
}

bool PolygonMesh::isSingularVertex(const int iV) const {
  int nV = getNumberOfVertices();
  return (0<=iV && iV<nV && _nPartsVertex[iV]>1);
}

// properties of the whole mesh

bool PolygonMesh::isRegular() const {
    int iV = 0; int iE = 0;
    int nE = getNumberOfEdges();
    int nV = getNumberOfVertices();
    while(isRegularEdge(iE++) && iE < nE );
    if(iE < nE) return false;
    while(!_isBoundaryVertex[iV++] && _nPartsVertex[iV] <=1 && iV < nV );
    if(iV < nV) return false;
    return true;
}

bool PolygonMesh::hasBoundary() const {
    int n = 0; int iE = 0; int nE = getNumberOfEdges();
    while((n+= isBoundaryEdge(iE++)) == 0 && iE < nE);
    return n > 0;
}
