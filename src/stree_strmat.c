#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "stree_strmat.h"

#define OPT_NODE_SIZE 24
#define OPT_LEAF_SIZE 12
#define OPT_INTLEAF_SIZE 12

SUFFIX_TREE stree_new_tree(int alphasize, int copyflag,
                           int build_type)
{
  SUFFIX_TREE tree;

  if (alphasize <= 0 || alphasize > 128)
    return NULL;

  if (build_type != LINKED_LIST &&
      build_type != SORTED_LIST &&
      build_type != COMPLETE_ARRAY)
    return NULL;

  if ((tree = malloc(sizeof(STREE_STRUCT))) == NULL)
    return NULL;

  memset(tree, 0, sizeof(STREE_STRUCT));
  tree->copyflag = copyflag;
  tree->alpha_size = alphasize;
  tree->build_type = build_type;

  if ((tree->root = int_stree_new_node(tree, NULL, NULL, 0)) == NULL) {
    free(tree);
    return NULL;
  }
  tree->num_nodes = 1;

  return tree;
}

void stree_delete_tree(SUFFIX_TREE tree)
{
  int i;

  int_stree_delete_subtree(tree, stree_get_root(tree));

  if (tree->strings != NULL) {
    if (tree->copyflag) {
      for (i=0; i < tree->strsize; i++)
        if (tree->strings[i] != NULL)
          free(tree->strings[i]);
    }
    free(tree->strings);
  }
  if (tree->rawstrings != NULL) {
    if (tree->copyflag) {
      for (i=0; i < tree->strsize; i++)
        if (tree->rawstrings[i] != NULL)
          free(tree->rawstrings[i]);
    }
    free(tree->rawstrings);
  }
  if (tree->ids != NULL)
    free(tree->ids);
  if (tree->lengths != NULL)
    free(tree->lengths);

  free(tree);
}

void stree_traverse(SUFFIX_TREE tree, 
                    int (*preorder_fn)(SUFFIX_TREE, STREE_NODE),
                    int (*postorder_fn)(SUFFIX_TREE, STREE_NODE))
{
  stree_traverse_subtree(tree, stree_get_root(tree), preorder_fn,
                         postorder_fn);
}

void stree_traverse_subtree(SUFFIX_TREE tree, STREE_NODE node,
                            int (*preorder_fn)(SUFFIX_TREE, STREE_NODE), 
                            int (*postorder_fn)(SUFFIX_TREE, STREE_NODE))
{
  enum { START, FIRST, MIDDLE, DONE, DONELEAF } state;
  int i, num, childnum;
  STREE_NODE root, child, *children;

  root = node;
  state = START;
  while (1) {
    if (state == START) {
      if (preorder_fn != NULL)
      (*preorder_fn)(tree, node);
      
      num = stree_get_num_children(tree, node);
      if (num > 0)
      state = FIRST;
      else
      state = DONELEAF;
      }
    
    if (state == FIRST || state == MIDDLE) {
      if (state == FIRST)
      childnum = 0;
      else
      childnum = node->isaleaf;
      
      if (!node->isanarray) {
        child = node->children;
        for (i=0; child != NULL && i < childnum; i++)
        child = child->next;
        }
      else {
        children = (STREE_NODE *) node->children;
        for (i=childnum; i < tree->alpha_size; i++) {
          if (children[i] != NULL)
          break;
          }
        child = (i < tree->alpha_size ? children[i] : NULL);
        }
      
      if (child == NULL)
      state = DONE;
      else {
        node->isaleaf = i + 1;
        node = child;
        
        state = START;
        }
      }
    
    if (state == DONE || state == DONELEAF) {
      if (state == DONE)
      node->isaleaf = 0;
      
      if (postorder_fn != NULL)
      (*postorder_fn)(tree, node);
      
      if (node == root)
      break;
      
      node = stree_get_parent(tree, node);
      state = MIDDLE;
      }
    
  }
}

int stree_match(SUFFIX_TREE tree, char *T, int N,
  STREE_NODE *node_out, int *pos_out)
{
  return stree_walk(tree, stree_get_root(tree), 0, T, N, node_out, pos_out);
  }

int stree_walk(SUFFIX_TREE tree, STREE_NODE node, int pos, char *T, int N,
STREE_NODE *node_out, int *pos_out)
{
  int len, endpos, edgelen;
  char *edgestr;
  STREE_NODE endnode;
  
  len = int_stree_walk_to_leaf(tree, node, pos, T, N, &endnode, &endpos);
  
  if (!int_stree_isaleaf(tree, endnode) || len == N) {
    *node_out = endnode;
    *pos_out = endpos;
    return len;
    }
  
  edgestr = stree_get_edgestr(tree, endnode);
  edgelen = stree_get_edgelen(tree, endnode);
  
  while (len < N && endpos < edgelen && T[len] == edgestr[endpos]) {
    len++;
    endpos++;
    }

  *node_out = endnode;
  *pos_out = endpos;
  return len;
}

STREE_NODE stree_find_child(SUFFIX_TREE tree, STREE_NODE node, char ch)
{
  char childch;
  STREE_NODE child, *children;

  if (int_stree_isaleaf(tree, node) || node->children == NULL)
    return NULL;

  if (!node->isanarray) {
    for (child=node->children; child != NULL; child=child->next) {
      childch = stree_getch(tree, child);

      if (ch == childch)
        return child;
      else if (tree->build_type == SORTED_LIST && ch < childch)
        return NULL;
    }

    return NULL;
  }
  else {
    children = (STREE_NODE *) node->children;

    return children[(int) ch];
  }
}

int stree_get_num_children(SUFFIX_TREE tree, STREE_NODE node)
{
  int i, count;
  STREE_NODE child, *children;

  if (int_stree_isaleaf(tree, node) || node->children == NULL)
    return 0;

  if (!node->isanarray) {
    count = 0;
    for (child=node->children; child != NULL; child=child->next)
      count++;
  }
  else {
    count = 0;
    children = (STREE_NODE *) node->children;
    for (i=0; i < tree->alpha_size; i++)
      if (children[i] != NULL)
        count++;
  }

  return count;
}

STREE_NODE stree_get_children(SUFFIX_TREE tree, STREE_NODE node)
{
  int i;
  STREE_NODE head, tail, *children;
  
  if (int_stree_isaleaf(tree, node) || node->children == NULL)
    return NULL;
  else if (!node->isanarray)
    return node->children;

  head = tail = NULL;
  children = (STREE_NODE *) node->children;
  for (i=0; i < tree->alpha_size; i++) {
    if (children[i] != NULL) {
      if (head == NULL)
        head = tail = children[i];
      else
        tail = tail->next = children[i];
    }
  }
  tail->next = NULL;

  return head;
}

int stree_get_num_leaves(SUFFIX_TREE tree, STREE_NODE node)
{
  int i;
  STREE_INTLEAF intleaf;

  if (int_stree_isaleaf(tree, node))
    return 1;
  else {
    intleaf = int_stree_get_intleaves(tree, node);
    for (i=0; intleaf != NULL; i++)
      intleaf = intleaf->next;
    return i;
  }
}

int stree_get_leaf(SUFFIX_TREE tree, STREE_NODE node, int leafnum,
                   char **string_out, int *pos_out, int *id_out)
{
  int i;
  STREE_INTLEAF intleaf;
  STREE_LEAF leaf;

  if (int_stree_isaleaf(tree, node)) {
    if (leafnum != 1)
      return 0;

    leaf = (STREE_LEAF) node;
    *string_out = int_stree_get_string(tree, leaf->strid);
    *pos_out = leaf->pos;
    *id_out = int_stree_get_strid(tree, leaf->strid);
    return 1;
  }
  else {
    intleaf = int_stree_get_intleaves(tree, node);
    for (i=1; intleaf != NULL; i++,intleaf=intleaf->next) {
      if (i == leafnum) {
        *string_out = int_stree_get_string(tree, intleaf->strid);
        *pos_out = intleaf->pos;
        *id_out = int_stree_get_strid(tree, intleaf->strid);
        return 1;
      }
    }
    return 0;
  }
}

void stree_reset_stats(SUFFIX_TREE tree)
{
  tree->num_compares = tree->edges_traversed = tree->links_traversed = 0;
  tree->child_cost = tree->nodes_created = tree->creation_cost = 0;
}

int stree_ukkonen_add_string(SUFFIX_TREE tree, char *S, char *Sraw,
                             int M, int strid)
{
  int i, j, g, h, gprime, edgelen, id;
  char *edgestr;
  STREE_NODE node, lastnode, root, child, parent;
  STREE_LEAF leaf;

  id = int_stree_insert_string(tree, S, Sraw, M, strid);
  if (id == -1)
    return 0;

  root = stree_get_root(tree);
  node = lastnode = root;
  g = 0;
  edgelen = 0;
  edgestr = NULL;

  for (i=0,j=0; i <= M; i++)  {
    for ( ; j <= i && j < M; j++) {
      if (g == 0 || g == edgelen) {
        if (i < M) {
          if ((child = stree_find_child(tree, node, S[i])) != NULL) {
            node = child;
            g = 1;
            edgestr = stree_get_edgestr(tree, node);
            edgelen = stree_get_edgelen(tree, node);
            break;
            }
          
          if ((leaf = int_stree_new_leaf(tree, id, i, j)) == NULL ||
          (node = int_stree_connect(tree, node,
          (STREE_NODE) leaf)) == NULL) {
            if (leaf != NULL)
            int_stree_free_leaf(tree, leaf);
            return 0;
            }
          
          tree->num_nodes++;
          }
        else {
          if (int_stree_isaleaf(tree, node) &&
          (node = int_stree_convert_leafnode(tree, node)) == NULL)
          return 0;
          
          if (!int_stree_add_intleaf(tree, node, id, j))
          return 0;
          }
        

        if (lastnode != root && lastnode->suffix_link == NULL)
          lastnode->suffix_link = node;
        lastnode = node;
      }
      else {
        if (i < M && S[i] == edgestr[g]) {
          g++;
          break;
        }

        if ((node = int_stree_edge_split(tree, node, g)) == NULL)
          return 0;

        edgestr = stree_get_edgestr(tree, node);
        edgelen = stree_get_edgelen(tree, node);

        if (i < M) {
          if ((leaf = int_stree_new_leaf(tree, id, i, j)) == NULL ||
          (node = int_stree_connect(tree, node,
          (STREE_NODE) leaf)) == NULL) {
            if (leaf != NULL)
            int_stree_free_leaf(tree, leaf);
            return 0;
            }
          
          tree->num_nodes++;
          }
        else {
          if (int_stree_isaleaf(tree, node) &&
          (node = int_stree_convert_leafnode(tree, node)) == NULL)
          return 0;
          
          if (!int_stree_add_intleaf(tree, node, id, j))
          return 0;
          }
        

        if (lastnode != root && lastnode->suffix_link == NULL)
          lastnode->suffix_link = node;
        lastnode = node;
      }

      if (node == root)
        ;
        else if (g == edgelen && node->suffix_link != NULL) {
          node = node->suffix_link;
          
          edgestr = stree_get_edgestr(tree, node);
          edgelen = stree_get_edgelen(tree, node);
          
          g = edgelen;
          continue;
          }
        else {
          parent = stree_get_parent(tree, node);
          if (parent != root) {
            node = parent->suffix_link;
            }
          else {
            node = root;
            g--;
            }
          edgelen = stree_get_edgelen(tree, node);
          
          h = i - g;
          while (g > 0) {
            node = stree_find_child(tree, node, S[h]);
            
            gprime = stree_get_edgelen(tree, node);
            if (gprime > g)
            break;
            
            g -= gprime;
            h += gprime;
            }
          
          edgelen = stree_get_edgelen(tree, node);
          edgestr = stree_get_edgestr(tree, node);
          
          if (g == 0) {
            if (lastnode != root && !int_stree_isaleaf(tree, node) &&
            lastnode->suffix_link == NULL) {
              lastnode->suffix_link = node;
              lastnode = node;
              }
      
            if (node != root)
            g = edgelen;
            }
          }
      
    }
  }

  return 1;
}

int int_stree_insert_string(SUFFIX_TREE tree, char *S, char *Sraw,
                            int M, int strid)
{
  int i, slot, newsize;

  if (tree->nextslot == tree->strsize) {
    if (tree->strsize == 0) {
      tree->strsize = 128;
      if ((tree->strings = malloc(tree->strsize * sizeof(char *))) == NULL ||
          (tree->rawstrings = malloc(tree->strsize * sizeof(char *))) == NULL)
        return -1;
      if ((tree->lengths = malloc(tree->strsize * sizeof(int))) == NULL ||
          (tree->ids = malloc(tree->strsize * sizeof(int))) == NULL)
        return -1;
      for (i=0; i < 128; i++) {
        tree->strings[i] = tree->rawstrings[i] = NULL;
        tree->lengths[i] = tree->ids[i] = 0;
      }
    }
    else {
      newsize = tree->strsize + tree->strsize;
      if ((tree->strings = realloc(tree->strings,
                                   tree->strsize * sizeof(char *))) == NULL ||
          (tree->rawstrings = realloc(tree->rawstrings,
                                      tree->strsize * sizeof(char *))) == NULL)
        return -1;
      if ((tree->lengths = realloc(tree->lengths,
                                   tree->strsize * sizeof(int))) == NULL ||
          (tree->ids = realloc(tree->ids,
                               tree->strsize * sizeof(int))) == NULL)
        return -1;

      for (i=tree->strsize; i < newsize; i++) {
        tree->strings[i] = tree->rawstrings[i] = NULL;
        tree->lengths[i] = tree->ids[i] = 0;
      }
      tree->strsize = newsize;
    }
  }

  slot = tree->nextslot;
  tree->strings[slot] = S;
  tree->rawstrings[slot] = Sraw;
  tree->lengths[slot] = M;
  tree->ids[slot] = strid;

  for (i=slot+1; i < tree->strsize; i++)
    if (tree->strings[i] == NULL)
      break;
  tree->nextslot = i;

  return slot;
}

void int_stree_delete_string(SUFFIX_TREE tree, int id)
{
  if (tree->strings[id] == NULL)
    return;

  if (tree->copyflag)
    free(tree->strings[id]);

  tree->strings[id] = NULL;
  if (id < tree->nextslot)
    tree->nextslot = id;
}

STREE_NODE int_stree_convert_leafnode(SUFFIX_TREE tree, STREE_NODE node)
{
  STREE_NODE newnode;
  STREE_LEAF leaf;
  STREE_INTLEAF ileaf;

  leaf = (STREE_LEAF) node;

  newnode = int_stree_new_node(tree, leaf->edgestr, leaf->rawedgestr,
                               leaf->edgelen);
  if (newnode == NULL)
    return NULL;

  if ((ileaf = int_stree_new_intleaf(tree, leaf->strid, leaf->pos)) == NULL) {
    int_stree_free_node(tree, newnode);
    return NULL;
  }

  newnode->id = leaf->id;
  newnode->leaves = ileaf;

  int_stree_reconnect(tree, node->parent, node, newnode);
  int_stree_free_leaf(tree, leaf);

  return newnode;
}

STREE_NODE int_stree_get_suffix_link(SUFFIX_TREE tree, STREE_NODE node)
{
  int len, edgelen;
  char *edgestr;
  STREE_NODE parent;

  if (node == stree_get_root(tree))
    return NULL;
  else if (!int_stree_isaleaf(tree, node))
    return node->suffix_link;

  edgestr = stree_get_edgestr(tree, node);
  edgelen = stree_get_edgelen(tree, node);
  parent = stree_get_parent(tree, node);

  if (parent != stree_get_root(tree))
    parent = parent->suffix_link;
  else {
    edgestr++;
    edgelen--;
  }

  node = parent;
  while (edgelen > 0) {
    node = stree_find_child(tree, node, *edgestr);
    assert(node != NULL);

    len = stree_get_edgelen(tree, node);
    edgestr += len;
    edgelen -= len;
  }

  return node;
}

STREE_NODE int_stree_connect(SUFFIX_TREE tree, STREE_NODE parent,
                             STREE_NODE child)
{
  int count;
  char ch;
  STREE_NODE temp, back, *children;

  if (int_stree_isaleaf(tree, parent) &&
      (parent = int_stree_convert_leafnode(tree, parent)) == NULL)
    return NULL;

  child->parent = parent;
  ch = stree_getch(tree, child);

  switch (tree->build_type) {
  case LINKED_LIST:
    child->next = parent->children;
    parent->children = child;
    break;

  case SORTED_LIST:
    back = NULL;
    for (temp=parent->children; temp != NULL; back=temp,temp=temp->next) {
      if (ch < stree_getch(tree, temp))
        break;
    }

    child->next = temp;
    if (back == NULL)
      parent->children = child;
    else
      back->next = child;
    break;

  case COMPLETE_ARRAY:
    children = (STREE_NODE *) parent->children;
    children[(int) ch] = child;
    break;
  }

  tree->idents_dirty = 1;

  return parent;
}

void int_stree_reconnect(SUFFIX_TREE tree, STREE_NODE parent,
                         STREE_NODE oldchild, STREE_NODE newchild)
{
  STREE_NODE child, back, *children;

  if (!parent->isanarray) {
    back = NULL;
    for (child=parent->children; child != oldchild; child=child->next)
      back = child;

    newchild->next = child->next;
    if (back == NULL)
      parent->children = newchild;
    else
      back->next = newchild;
  }
  else {
    children = (STREE_NODE *) parent->children;
    children[(int) stree_getch(tree, newchild)] = newchild;
  }

  newchild->parent = parent;
  oldchild->parent = NULL;

  tree->idents_dirty = 1;
}

void int_stree_disc_from_parent(SUFFIX_TREE tree, STREE_NODE parent,
                                STREE_NODE child)
{
  STREE_NODE node, back, *children;

  if (!parent->isanarray) {
    back = NULL;
    for (node=parent->children; node != child; node=node->next)
      back = node;

    if (back == NULL)
      parent->children = node->next;
    else
      back->next = node->next;
  }
  else {
    children = (STREE_NODE *) parent->children;
    children[(int) stree_getch(tree, child)] = NULL;
  }
}

void int_stree_disconnect(SUFFIX_TREE tree, STREE_NODE node)
{
  int num;
  STREE_NODE parent;

  if (node == stree_get_root(tree))
    return;

  parent = stree_get_parent(tree, node);
  int_stree_disc_from_parent(tree, parent, node);

  if (parent->leaves == NULL && parent != stree_get_root(tree) && 
      (num = stree_get_num_children(tree, parent)) < 2) {
    if (num == 0) {
      int_stree_disconnect(tree, parent);
      int_stree_delete_subtree(tree, parent);
    }
    else if (num == 1)
      int_stree_edge_merge(tree, parent);
  }

  tree->idents_dirty = 1;
}

STREE_NODE int_stree_edge_split(SUFFIX_TREE tree, STREE_NODE node, int len)
{
  STREE_NODE newnode, parent;

  if (node == stree_get_root(tree) ||
      len == 0 || stree_get_edgelen(tree, node) <= len)
    return NULL;

  newnode = int_stree_new_node(tree, node->edgestr, node->rawedgestr, len);
  if (newnode == NULL)
    return NULL;

  parent = stree_get_parent(tree, node);
  int_stree_reconnect(tree, parent, node, newnode);
  
  node->edgestr += len;
  node->rawedgestr += len;
  node->edgelen -= len;

  if (int_stree_connect(tree, newnode, node) == NULL) {
    node->edgestr -= len;
    node->rawedgestr -= len;
    node->edgelen += len;
    int_stree_reconnect(tree, parent, newnode, node);
    int_stree_free_node(tree, newnode);
    return NULL;
  }

  tree->num_nodes++;
  tree->idents_dirty = 1;

  return newnode;
}

void int_stree_edge_merge(SUFFIX_TREE tree, STREE_NODE node)
{
  int i, len;
  STREE_NODE parent, child, *children;

  if (node == stree_get_root(tree) || int_stree_isaleaf(tree, node) ||
      node->leaves != NULL)
    return;
  
  parent = stree_get_parent(tree, node);
  if (!node->isanarray) {
    child = node->children;
    if (child == NULL || child->next != NULL)
      return;
  }
  else {
    child = NULL;
    children = (STREE_NODE *) node->children;
    for (i=0; i < tree->alpha_size; i++) {
      if (children[i] != NULL) {
        if (child != NULL)
          return;
        child = children[i];
      }
    }
    if (child == NULL)
      return;
  }
  len = stree_get_edgelen(tree, node);
  child->edgestr -= len;
  child->rawedgestr -= len;
  child->edgelen += len;

  int_stree_reconnect(tree, parent, node, child);
  tree->num_nodes--;
  tree->idents_dirty = 1;

  int_stree_free_node(tree, node);
}

int int_stree_add_intleaf(SUFFIX_TREE tree, STREE_NODE node,
                          int strid, int pos)
{
  STREE_INTLEAF intleaf;

  if (int_stree_isaleaf(tree, node) ||
      (intleaf = int_stree_new_intleaf(tree, strid, pos)) == NULL)
    return 0;

  intleaf->next = node->leaves;
  node->leaves = intleaf;
  return 1;
}

int int_stree_remove_intleaf(SUFFIX_TREE tree, STREE_NODE node,
                             int strid, int pos)
{
  STREE_INTLEAF intleaf, back;

  if (int_stree_isaleaf(tree, node) || !int_stree_has_intleaves(tree, node))
    return 0;

  back = NULL;
  intleaf = int_stree_get_intleaves(tree, node);
  for (back=NULL; intleaf != NULL; back=intleaf,intleaf=intleaf->next)
    if (intleaf->strid == strid && intleaf->pos == pos)
      break;

  if (intleaf == NULL)
    return 0;

  if (back != NULL)
    back->next = intleaf->next;
  else
    node->leaves = intleaf->next;

  int_stree_free_intleaf(tree, intleaf);
  return 1;
}

void int_stree_delete_subtree(SUFFIX_TREE tree, STREE_NODE node)
{
  int i;
  STREE_NODE child, temp, *children;
  STREE_INTLEAF ileaf, itemp;

  if (int_stree_isaleaf(tree, node))
    int_stree_free_leaf(tree, (STREE_LEAF) node);
  else {
    for (ileaf=node->leaves; ileaf != NULL; ileaf=itemp) {
      itemp = ileaf->next;
      int_stree_free_intleaf(tree, ileaf);
    }

    if (!node->isanarray) {
      for (child=node->children; child != NULL; child=temp) {
        temp = child->next;
        int_stree_delete_subtree(tree, child);
      }
    }
    else {
      children = (STREE_NODE *) node->children;
      for (i=0; i < tree->alpha_size; i++)
        if (children[i] != NULL)
          int_stree_delete_subtree(tree, children[i]);
    }

    int_stree_free_node(tree, node);
  }
}

int int_stree_walk_to_leaf(SUFFIX_TREE tree, STREE_NODE node, int pos,
                           char *T, int N, STREE_NODE *node_out, int *pos_out)
{
  int len, edgelen;
  char *edgestr;
  STREE_NODE child;

  if (int_stree_isaleaf(tree, node)) {
    *node_out = node;
    *pos_out = pos;
    return 0;
  }

  edgestr = stree_get_edgestr(tree, node);
  edgelen = stree_get_edgelen(tree, node);
  len = 0;
  while (1) {
    while (len < N && pos < edgelen && T[len] == edgestr[pos]) {
      pos++;
      len++;
    }

    if (len == N || pos < edgelen ||
        (child = stree_find_child(tree, node, T[len])) == NULL)
      break;

    if (int_stree_isaleaf(tree, child)) {
      *node_out = child;
      *pos_out = 0;
      return len;
    }

    node = child;
    edgestr = stree_get_edgestr(tree, node);
    edgelen = stree_get_edgelen(tree, node);
    pos = 1;
    len++;
  }

  *node_out = node;
  *pos_out = pos;
  return len;
}

void int_stree_set_idents(SUFFIX_TREE tree)
{
  enum { START, FIRST, MIDDLE, DONE, DONELEAF } state;
  int i, num, childnum, nextid;
  STREE_NODE root, node, child, *children;

  if (!tree->idents_dirty)
    return;

  nextid = 0;
  node = root = stree_get_root(tree);
  state = START;
  while (1) {
    if (state == START) {
        node->id = nextid++;
        
        num = stree_get_num_children(tree, node);
        if (num > 0)
        state = FIRST;
        else
        state = DONELEAF;
        }
    
    if (state == FIRST || state == MIDDLE) {
        if (state == FIRST)
        childnum = 0;
        else
        childnum = node->isaleaf;
        
        if (!node->isanarray) {
          child = node->children;
          for (i=0; child != NULL && i < childnum; i++)
          child = child->next;
          }
        else {
          children = (STREE_NODE *) node->children;
          for (i=childnum; i < tree->alpha_size; i++)
          if (children[i] != NULL)
          break;
          child = (i < tree->alpha_size ? children[i] : NULL);
          }
        
        if (child == NULL)
        state = DONE;
        else {
          node->isaleaf = i + 1;
          node = child;
          state = START;
          }
        }
    
    if (state == DONE || state == DONELEAF) {
        if (state == DONE)
        node->isaleaf = 0;
        
        if (node == root)
        break;
        
        node = node->parent;
        state = MIDDLE;
        }
    
  }

  tree->idents_dirty = 0;
}

STREE_INTLEAF int_stree_new_intleaf(SUFFIX_TREE tree, int strid, int pos)
{
  STREE_INTLEAF ileaf;
  
  if ((ileaf = malloc(sizeof(SINTLEAF_STRUCT))) == NULL)
    return NULL;

  memset(ileaf, 0, sizeof(SINTLEAF_STRUCT));
  ileaf->strid = strid;
  ileaf->pos = pos;

#ifdef STATS
  tree->tree_size += OPT_INTLEAF_SIZE;
#endif

  return ileaf;
}

STREE_LEAF int_stree_new_leaf(SUFFIX_TREE tree, int strid, int edgepos,
  int leafpos)
{
  STREE_LEAF leaf;

  if ((leaf = malloc(sizeof(SLEAF_STRUCT))) == NULL)
    return NULL;

  memset(leaf, 0, sizeof(SLEAF_STRUCT));
  leaf->isaleaf = 1;
  leaf->strid = strid;
  leaf->pos = leafpos;
  leaf->edgestr = int_stree_get_string(tree, strid) + edgepos;
  leaf->rawedgestr = int_stree_get_rawstring(tree, strid) + edgepos;
  leaf->edgelen = int_stree_get_length(tree, strid) - edgepos;

#ifdef STATS
  tree->tree_size += OPT_LEAF_SIZE;
#endif

  return leaf;
}

STREE_NODE int_stree_new_node(SUFFIX_TREE tree, char *edgestr,
  char *rawedgestr, int edgelen)
{
  STREE_NODE node;

  if ((node = malloc(sizeof(SNODE_STRUCT))) == NULL)
    return NULL;

  memset(node, 0, sizeof(SNODE_STRUCT));
  node->edgestr = edgestr;
  node->rawedgestr = rawedgestr;
  node->edgelen = edgelen;

  if (tree->build_type == COMPLETE_ARRAY) {
    node->children = malloc(tree->alpha_size * sizeof(STREE_NODE));
    if (node->children == NULL) {
      free(node);
      return NULL;
    }

    memset(node->children, 0, tree->alpha_size * sizeof(STREE_NODE));
    node->isanarray = 1;

#ifdef STATS
    tree->tree_size += tree->alpha_size * sizeof(STREE_NODE);
#endif
  }

#ifdef STATS
  tree->tree_size += OPT_NODE_SIZE;
#endif

  return node;
}

void int_stree_free_intleaf(SUFFIX_TREE tree, STREE_INTLEAF ileaf)
{
#ifdef STATS
  tree->tree_size -= OPT_INTLEAF_SIZE;
#endif

  free(ileaf);
}

void int_stree_free_leaf(SUFFIX_TREE tree, STREE_LEAF leaf)
{
#ifdef STATS
  tree->tree_size -= OPT_LEAF_SIZE;
#endif

  free(leaf);
}

void int_stree_free_node(SUFFIX_TREE tree, STREE_NODE node)
{
  if (node->isanarray) {
    free(node->children);

#ifdef STATS
    tree->tree_size -= tree->alpha_size * sizeof(STREE_NODE);
#endif
  }

#ifdef STATS
  tree->tree_size -= OPT_NODE_SIZE;
#endif
  free(node);
}

