#ifndef MARK_H
#define MARK_H

#include "kdNode.h"
#include "pbbs/unionFind.h"

template<class nodeT>
void markAll(nodeT *nd, intT id){
  if(!nd->isLeaf() && nd->getId() != id){
    if (nd->size() > 2000) {
      par_do([&](){markAll(nd->L(), id);},
	     [&](){markAll(nd->R(), id);});
    } else {
      markAll(nd->L(), id);
      markAll(nd->R(), id);}
  }
  nd->setId(id);
}

template<class nodeT, class pointT>
void mark(nodeT *nd, edgeUnionFind *uf, pointT* s){//todo look for that optimization

  if(nd->hasId()) {
    return markAll(nd, uf->find(nd->getItems()[0]-s));
  }

  nd->setId(uf->find(nd->getItems()[0]-s));

  if(nd->isLeaf()){
    for (intT i=1; i < nd->size(); ++i) {
      if (nd->getId() != uf->find(nd->getItems()[i]-s)) {
        return nd->resetId();}
    }
  } else {
    if (nd->size() > 2000) {
      par_do([&](){mark(nd->L(), uf, s);},
	     [&](){mark(nd->R(), uf, s);});
    } else {
      mark(nd->L(), uf, s);
      mark(nd->R(), uf, s);
    }
    if (nd->getId() != nd->L()->getId()) return nd->resetId();
    if (nd->getId() != nd->R()->getId()) return nd->resetId();
  }
}

#endif
