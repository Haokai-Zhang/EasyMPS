#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Author: Hao-Kai Zhang <zhk20@tsinghua.mails.edu.cn>
Creation Date: 2020-07-01 11:00
Description: EasyMPS project. <class_fsa.py> contains finite-state automata generating matrix product operator automatically.
'''

import numpy as np
from AutoMPO.mpo_gadgets import IsClose, IsOrdered, BoundSort
from AutoMPO.class_named_data import named_data

class fsa_edge(object):

    def __init__(self, named_opr, weight, node_from, node_to):
        '''
        initialization of edge class of a finite state automata
        :param named_opr: named operator on this edge
        :param weight: weight on this edge
        :param node_from: start node of the edge
        :param node_to: end node of the edge
        '''
        # the attachments on each edge: {operator, weight}
        self.named_opr = named_opr
        self.weight = weight
        # the edge from node_from to node_to
        # (node_from) -----edge----> (node_to)
        self.node_from = node_from
        self.node_to = node_to
        # this edge is the out_idx-th out-edge of node_from
        self.out_idx = node_from.OutDegree()
        # this edge is the in_idx-th in-edge of node_to
        self.in_idx = node_to.InDegree()

    def EqualEdge(self, other):
        # two edge is equal to each other if they have the same attachments
        # here we do not use "__eq__" for later convenience
        return (self.named_opr == other.named_opr) and IsClose(self.weight, other.weight)

    def IsWanted(self, wanted_named_opr, wanted_weight):
        # decide whether the edge is wanted by its attachments
        return (self.named_opr == wanted_named_opr) and IsClose(self.weight, wanted_weight)


class fsa_node(object):

    def __init__(self, layer_idx, bond_idx):
        '''
        initialization of node class of a finite state automata
        :param layer_idx: space index of the node
        :param bond_idx: bond index of the node
        '''
        # layer_idx and bond_idx record the location of this node on the whole graph
        # layer_idx records which layer it belongs to
        self.layer_idx = layer_idx
        # bond_idx records the serial number in the layer
        self.bond_idx = bond_idx
        # create two lists to record the out-edges and in-edges from/to this node
        self.list_edge_out = []
        self.list_edge_in = []

    def __eq__(self, other):
        '''
        decide whether two nodes are identical
        :param other: the other node
        :return: whether two nodes are identical
        '''
        # two nodes are identical if they have the same location
        return (self.layer_idx == other.layer_idx) and (self.bond_idx == other.bond_idx)

    def OutDegree(self):
        '''
        :return: the out-degree of this node
        '''
        return len(self.list_edge_out)

    def InDegree(self):
        '''
        :return: the in-degree of this node
        '''
        return len(self.list_edge_in)

    def NewEdgeOut(self, next_node, named_opr, weight):
        '''
        create a new out-edge to an existed node
        :param next_node: the existed next node
        :param named_opr: the operator on this new edge
        :param weight: the weight on this new edge
        :return: the new out-edge
        '''
        # create a new edge from self to next_node
        new_edge = fsa_edge(named_opr, weight, self, next_node)
        # record this edge on the list_edge_out of self
        self.list_edge_out.append(new_edge)
        # record this edge on the list_edge_in of next_node
        next_node.list_edge_in.append(new_edge)
        return new_edge

    def SearchEdgeOut(self, wanted_named_opr, wanted_weight):
        '''
        search for the out-edge from this node
        :param wanted_named_opr: wanted operator
        :param wanted_weight: wanted weight
        :return: the first out-edge that satisfies the wanted condition
        '''
        for the_edge in self.list_edge_out:
            if the_edge.IsWanted(wanted_named_opr, wanted_weight):
                return the_edge
        return None

    def SearchEdgeIn(self, wanted_named_opr, wanted_weight):
        '''
        search for the in-edge from this node
        :param wanted_named_opr: wanted operator
        :param wanted_weight: wanted weight
        :return: the first in-edge that satisfies the wanted condition
        '''
        for the_edge in self.list_edge_in:
            if the_edge.IsWanted(wanted_named_opr, wanted_weight):
                return the_edge
        return None

    def SearchHostEdgeOut(self, merged_edge):
        '''
        search for the host (out-)edge to merge into
        :param merged_edge: the (out-)edge to be merged
        :return: the host edge, if any
        '''
        for the_edge in self.list_edge_out:
            if merged_edge.EqualEdge(the_edge) and (the_edge.node_to != merged_edge.node_to):
                return the_edge
        return None

    def SearchHostEdgeIn(self, merged_edge):
        '''
        search for the host (in-)edge to merge into
        :param merged_edge: the (in-)edge to be merged
        :return: the host edge, if any
        '''
        for the_edge in self.list_edge_in:
            if merged_edge.EqualEdge(the_edge) and (the_edge.node_from != merged_edge.node_from):
                return the_edge
        return None

    def SearchOprOut(self, wanted_named_opr):
        '''
        search for the out-edge with the wanted operator from this node
        :param wanted_named_opr: wanted operator
        :return: the first out-edge with the wanted operator
        '''
        for the_edge in self.list_edge_out:
            if (the_edge.named_opr == wanted_named_opr):
                return the_edge
        return None


class fsa(object):

    def __init__(self, site_num, identity_name='Id'):

        # record the site number, matrix number = site number
        # node-layer number = site_number + 1, i.e.{0,1,2,...,site_num}
        self.site_num = site_num
        # create the head node and tail node
        self.head = fsa_node(layer_idx=0, bond_idx=0)
        self.tail = fsa_node(layer_idx=site_num, bond_idx=0)
        # create a 2-d list to point to each node to record the location
        self.locator = [[] for i in range(0, site_num + 1)]
        # record the location of head and tail
        self.locator[0].append(self.head)
        self.locator[-1].append(self.tail)
        # locator : ------------------------------------
        #           |                                  |
        #           |   /                           \  |
        #          head --         ......          --tail
        #               \                           /
        # record the name of the identity used in SortLayer()
        self.name_Id = identity_name

    def Add(self, coef, list_named_phys_opr, list_phys_site, list_named_inst_opr=None, print_form=None):
        '''
        add a new term to this finite state automata
        :param coef: coefficient corresponding to this term
        :param list_named_phys_opr: list of named physical operators
        :param list_phys_site: list of site serial number corresponding to the physical operators
        :param list_named_inst_opr: list of inserted operators
        :param print_form: the form of print = {None, 'phys', 'all'}
        :return: None
        '''
        # operator configuration:
        # inst_opr s --- phys_opr s--- inst_opr s -...- phys_opr s --- inst_opr s
        if (list_named_inst_opr == None):
            # default None means inserted operators are all identities
            opr_shape = list_named_phys_opr[0].data.shape
            Id = named_data(self.name_Id, np.identity(opr_shape[0]))
            list_named_inst_opr = [Id for i in range(0, len(list_named_phys_opr) + 1)]

        # assert whether the lengths of lists match
        assert len(list_named_phys_opr) == len(list_phys_site), '----Error: unmatched list length!----'
        assert len(list_named_phys_opr) + 1 == len(list_named_inst_opr), '----Error: unmatched list length!----'
        # assert all the phys_sites are smaller than the system size
        for phys_site in list_phys_site:
            assert phys_site < self.site_num, '----Error: site exceeding system size!----'

        # sort list_phys_site and list_named_phys_opr according to list_phys_site
        if not IsOrdered(list_phys_site):
            list_phys_site, list_named_phys_opr = BoundSort(list_phys_site, list_named_phys_opr)
            print('----Warning: physical operators are not in order, '
                  'which have been sorted automatically without changing inserted operators.----')

        # default weight = 1
        list_weight = [1. for i in range(0, self.site_num)]
        # put the coefficient on the weight of the first physical operator (convention)
        if (len(list_phys_site) == 0):
            list_weight[0] = coef
        else:
            list_weight[list_phys_site[0]] = coef

        # put the physical and inserted operators in order in list_named_opr
        list_named_opr = []
        # record the used physical operators
        idx_used_phys = 0
        for site in range(0, self.site_num):
            if (idx_used_phys < len(list_phys_site)) and (site == list_phys_site[idx_used_phys]):
                # if the physical operators are not used up
                # and there is a physical operator in the list corresponding to this site
                list_named_opr.append(list_named_phys_opr[idx_used_phys])
                idx_used_phys += 1
            else:
                # filled with inserted operators
                list_named_opr.append(list_named_inst_opr[idx_used_phys])

        # print the term just added
        if (print_form == 'all'):
            print(coef, end=' ')
            # print physical and inserted operators in order
            for site in range(0, self.site_num):
                print(list_named_opr[site].name, end=' ')
            print('')
        elif (print_form == 'phys'):
            print(coef, end=' ')
            # print physical operators with their site indices
            for phys_idx in range(0, len(list_named_phys_opr)):
                print(list_named_phys_opr[phys_idx].name + '_' + str(list_phys_site[phys_idx]), end=' ')
            print('')

        # create a new path according to the operator-weight list corresponding to each edge
        self.NewPath(list_named_opr, list_weight)

    def NewPath(self, list_named_opr, list_weight):
        '''
        create a new path for each term
        :param list_named_opr: list of named operators of the term
        :param list_weight: list of weights of each operator
        :return: None
        '''
        # assert whether the number of operators equals to the number of sites
        assert len(list_named_opr) == self.site_num, '----Error: illegal list length!----'
        # use a list to record the newly created edges for later convenience
        new_edge_register = []

        # combine like terms
        if self.CombineLikeTerm(list_named_opr, list_weight):
            return None

        # start from the head node
        last_node = self.head
        # node = 0,1,2,...,n-2,n-1
        # edge = 0,1,2,...,n-2
        for layer_idx in range(0, self.site_num - 1):
            # create a new node and a new edge
            the_node = self.NewSuccessor(last_node, list_named_opr[layer_idx], list_weight[layer_idx])
            # record this new edge in new_edge_register
            new_edge_register.append(the_node.list_edge_in[0])
            # proceed to the next node
            last_node = the_node

        # the tail node
        # node = n-1,n
        # edge = n-1
        layer_idx = self.site_num - 1
        # create a new edge to the tail node
        the_edge = self.NewEdge(last_node, self.tail, list_named_opr[layer_idx], list_weight[layer_idx])
        new_edge_register.append(the_edge)


        # forward merge
        # start from the head node to the tail node to check if there are edges that can be merged
        for the_new_edge in new_edge_register:
            if (the_new_edge.node_from.InDegree() > 1):
                # avoid cross-edge
                break
            else:
                # search for an old edge as a host edge to merge to
                host_edge = the_new_edge.node_from.SearchHostEdgeOut(the_new_edge)
                if (host_edge == None) or (host_edge.node_to.InDegree() > 1):
                    break
                else:
                    # merge the redundant edges and nodes
                    is_merged = self.MergeSuccessor(the_new_edge, host_edge)
                    if not is_merged:
                        break
        # backward merge
        # start from the tail node to the head node to check if there are edges that can be merged
        new_edge_register.reverse()
        for the_new_edge in new_edge_register:
            if (the_new_edge.node_from.OutDegree() > 1):
                # avoid cross-edge
                break
            else:
                # search for an old edge as a host edge to merge to
                host_edge = the_new_edge.node_to.SearchHostEdgeIn(the_new_edge)
                if (host_edge == None) or (host_edge.node_from.OutDegree() > 1):
                    break
                else:
                    # merge the redundant edges and nodes
                    is_merged = self.MergePrecursor(the_new_edge, host_edge)
                    if not is_merged:
                        break

    def CombineLikeTerm(self, list_named_opr, list_weight):
        '''
        combine like terms when we create a new path
        :param list_named_opr: list of named operators of the term
        :param list_weight: list of weights of each operator
        :return: True = there are like terms and they have been combined, False = there are no like terms
        '''
        # start from the head node
        last_node = self.head
        # record the weighted edge (convention: first physical edge)
        first_phys_edge = None
        for layer_idx in range(0, self.site_num):
            # search for the edge with the same operator
            the_edge = last_node.SearchOprOut(list_named_opr[layer_idx])
            if (the_edge == None):
                return False
            else:
                # if yes, continue to walk along this edge
                last_node = the_edge.node_to
                # check whether it is the weighted edge
                if (list_named_opr[layer_idx].name != self.name_Id) and (first_phys_edge == None):
                    first_phys_edge = the_edge
        # there are like terms
        # combine them by adding the weights together
        first_phys_edge.weight += list_weight[first_phys_edge.node_from.layer_idx]
        # raise warning to remind the like terms have been combined automatically
        print('----Warning: there are like terms, which have been combined automatically.----')
        return True

    def NewEdge(self, node_from, node_to, named_opr, weight):
        '''
        create a new edge
        :param node_from: new edge from node_from
        :param node_to: new edge to node_to
        :param named_opr: the operator on this new edge
        :param weight: the weight on this new edge
        :return: the created new edge
        '''
        new_edge = node_from.NewEdgeOut(node_to, named_opr, weight)
        return new_edge

    def NewSuccessor(self, last_node, named_opr, weight):
        '''
        create a new successor node
        :param last_node: the node which the successor node succeeds
        :param named_opr: the operator on the edge between last_node and successor node
        :param weight: the weight on the edge between last_node and successor node
        :return: the created new successor node
        '''
        # the successor locates at the next layer
        new_node_layer_idx = last_node.layer_idx + 1
        # the successor locates at the bottom of the next layer
        new_node_bond_idx = len(self.locator[new_node_layer_idx])
        # create a new node as the successor
        new_node = fsa_node(new_node_layer_idx, new_node_bond_idx)
        # create a new edge from the last_node to the new successor node
        self.NewEdge(last_node, new_node, named_opr, weight)
        # use locator to point to this new successor
        self.locator[new_node_layer_idx].append(new_node)
        return new_node

    def PrecursorCoincide(self, node0, node1):
        '''
        decide whether the precursors of two nodes coincide
        :param node0: one of the node
        :param node1: another one of the node
        :return: True=coincide, False=not coincide
        '''
        for edge0 in node0.list_edge_in:
            for edge1 in node1.list_edge_in:
                if (edge0.node_from == edge1.node_from):
                    return True
        return False

    def SuccessorCoincide(self, node0, node1):
        '''
        decide whether the successors of two nodes coincide
        :param node0: one of the node
        :param node1: another one of the node
        :return: True=coincide, False=not coincide
        '''
        for edge0 in node0.list_edge_out:
            for edge1 in node1.list_edge_out:
                if (edge0.node_to == edge1.node_to):
                    return True
        return False

    def MergePrecursor(self, merged_edge, host_edge):
        '''
        merge the precursors provided that they have the same successor
        :param merged_edge: the out-edge to be merged
        :param host_edge: the out-edge to merge to
        :return: whether the merge is successful
        '''
        # merging example
        # 0---|''|                        0---|''|
        # 1---|..|---.2(host)             1---|  |------.2(4)
        #             \        merge      3---|..|       \
        #              |''|   ------->                    |''|
        #              |..|                               |..|
        #             /
        # 3---|''|---'4(merged)
        #     |..|

        # the node_from of host/merged edge is host/merged node
        host_node = host_edge.node_from
        merged_node = merged_edge.node_from

        # assert the two nodes are in the same layer
        assert host_node.layer_idx == merged_node.layer_idx, '----Error: illegal layer_idx!----'
        # assert the out-degree of the host_node equals 1
        assert host_node.OutDegree() <= 1, '----Error: unexpected cross-edge!----'
        assert merged_node.OutDegree() <= 1, '----Error: unexpected cross-edge!----'
        # cross-edge example
        # "4-2"
        # 0---|''|---2                    0---|''|---2
        # 1---|..|---.3(host)             1---|  |------.3(5)
        #             \        merge      4---|..|       \
        #              |''|   ------->                    |''|
        #              |..|                               |..|
        #             /
        # 4---|''|---'5(merged)
        #     |..|

        if self.PrecursorCoincide(host_node, merged_node):
            # if the two nodes have the same precursor, they can not be merged
            # fail to merge
            return False
        else:
            # delete the merged_edge
            self.DelBackwardMergedEdge(merged_edge)
            # deliver edge_in from merged_node to host_node
            self.DeliverEdgeIn(merged_node, host_node)
            # change bond_idx larger than merged_node.bond_idx
            for the_node in self.locator[merged_node.layer_idx][(merged_node.bond_idx+1):]:
                the_node.bond_idx -= 1
            # delete the merged_node from the locator
            del self.locator[merged_node.layer_idx][merged_node.bond_idx]
            del merged_node
            # success to merge
            return True

    def MergeSuccessor(self, merged_edge, host_edge):
        '''
        merge the successor provided that they have the same precursor
        :param merged_edge: the in-edge to be merged
        :param host_edge: the in-edge to merge to
        :return: whether the merge is successful
        '''

        # the node_to of host/merged edge is host/merged node
        host_node = host_edge.node_to
        merged_node = merged_edge.node_to

        # assert the two nodes are in the same layer
        assert host_node.layer_idx == merged_node.layer_idx, '----Error: illegal layer_idx!----'
        # assert the in-degree of the host_node equals 1
        assert host_node.InDegree() <= 1, '----Error: unexpected cross-edge!----'
        assert merged_node.InDegree() <= 1, '----Error: unexpected cross-edge!----'

        if self.SuccessorCoincide(host_node, merged_node):
            # if the two nodes have the same successor, they can not be merged
            # fail to merge
            return False
        else:
            # delete the merged_edge
            self.DelForwardMergedEdge(merged_edge)
            # deliver edge_out from merged_node to host_node
            self.DeliverEdgeOut(merged_node, host_node)
            # change bond_idx larger than merged_node.bond_idx
            for the_node in self.locator[merged_node.layer_idx][(merged_node.bond_idx+1):]:
                the_node.bond_idx -= 1
            # delete the merged_node from the locator
            del self.locator[merged_node.layer_idx][merged_node.bond_idx]
            del merged_node
            # success to merge
            return True

    def DelForwardMergedEdge(self, merged_edge):
        '''
        delete the forward merged edge
        :param merged_edge: the edge to be merged
        :return: None
        '''
        # change edge's register with in_idx larger than merged_edge
        for other_edge_out in merged_edge.node_from.list_edge_out[(merged_edge.out_idx+1):]:
            other_edge_out.out_idx -= 1
        # delete merged_edge from merged_edge.node_to's register
        del merged_edge.node_from.list_edge_out[merged_edge.out_idx]
        # delete merged_node's register for its in-edges
        del merged_edge.node_to.list_edge_in[:]

    def DelBackwardMergedEdge(self, merged_edge):
        '''
        delete the backward merged edge
        :param merged_edge: the edge to be merged
        :return: None
        '''
        # change edge's register with in_idx larger than merged_edge
        for other_edge_in in merged_edge.node_to.list_edge_in[(merged_edge.in_idx+1):]:
            other_edge_in.in_idx -= 1
        # delete merged_edge from merged_edge.node_to's register
        del merged_edge.node_to.list_edge_in[merged_edge.in_idx]
        # delete merged_node's register for its out-edges
        del merged_edge.node_from.list_edge_out[:]

    def DeliverEdgeOut(self, merged_node, host_node):
        '''
        deliver all out-edges from merged_node to host_node
        :param merged_node: the node to be merged
        :param host_node: the node to merge to
        :param merged_edge: the edge to be merged
        :return: None
        '''
        for the_edge in merged_node.list_edge_out:
            # pull out all the out-edges from merged_node
            # plug into host_node
            # change edge's register to record node
            the_edge.node_from = host_node
            the_edge.out_idx = host_node.OutDegree()
            # change node's register to record edge
            host_node.list_edge_out.append(the_edge)
        # delete merged_node's register for its out-edges
        del merged_node.list_edge_out[:]

    def DeliverEdgeIn(self, merged_node, host_node):
        '''
        deliver all in-edges from merged_node to host_node
        :param merged_node: the node to be merged
        :param host_node: the node to merge to
        :param merged_edge: the edge to be merged
        :return: None
        '''
        for the_edge in merged_node.list_edge_in:
            # pull out all the in-edges from merged_node
            # plug into host_node
            # change edge's register to record node
            the_edge.node_to = host_node
            the_edge.in_idx = host_node.InDegree()
            # change node's register to record edge
            host_node.list_edge_in.append(the_edge)
        # delete merged_node's register for its in-edges
        del merged_node.list_edge_in[:]


    def GenMPO(self, sort_order=None):
        '''
        generate matrix product operators according to the finite state automata
        :param sort_order: the node with the identity edge from the head node is at the "bottom"/"top" of the layer
        :return: list of MPO
        '''
        # sort the order of nodes in each layer before generating
        if (sort_order != None):
            self.SortLayer(sort_order)
        # create a all-zero list to record MPO
        list_mpo = [[] for i in range(0, self.site_num)]
        opr_shape = self.head.list_edge_out[0].named_opr.data.shape
        for layer_idx in range(0, self.site_num):
            for i in range(0, len(self.locator[layer_idx])):
                list_mpo[layer_idx].append([])
                for j in range(0, len(self.locator[layer_idx + 1])):
                    list_mpo[layer_idx][i].append(np.zeros(opr_shape))
        # for each edge
        # row index i = node_from.bond_idx
        # column index j = node_to.bond_idx
        for layer_idx in range(0, self.site_num):
            for the_node in self.locator[layer_idx]:
                for the_edge in the_node.list_edge_out:
                    i = the_edge.node_from.bond_idx
                    j = the_edge.node_to.bond_idx
                    # matrix element = weight * local operator
                    list_mpo[layer_idx][i][j] = the_edge.weight * the_edge.named_opr.data
            list_mpo[layer_idx] = np.array(list_mpo[layer_idx])
        return list_mpo

    def GenSymbolMPO(self, sort_order='bottom'):
        '''
        generate matrix product operators at each site with local operators represented as just symbols
        :param sort_order: the node with the identity edge from the head node is at the "bottom"/"top" of the layer
        :return: list of symbol MPO
        '''
        # sort the order of nodes in each layer before generating
        if (sort_order != None):
            self.SortLayer(sort_order)
        # create a all-zero list to record the symbol MPO
        list_symbol_mpo = [[] for i in range(0, self.site_num)]
        for layer_idx in range(0, self.site_num):
            for i in range(0, len(self.locator[layer_idx])):
                list_symbol_mpo[layer_idx].append([])
                for j in range(0, len(self.locator[layer_idx + 1])):
                    list_symbol_mpo[layer_idx][i].append('0')
        # for each edge
        # row index i = node_from.bond_idx
        # column index j = node_to.bond_idx
        for layer_idx in range(0, self.site_num):
            for the_node in self.locator[layer_idx]:
                for the_edge in the_node.list_edge_out:
                    i = the_edge.node_from.bond_idx
                    j = the_edge.node_to.bond_idx
                    # matrix element = weight * local operator
                    list_symbol_mpo[layer_idx][i][j] = str(the_edge.weight) + the_edge.named_opr.name
        return list_symbol_mpo

    def PrintSymbolMPO(self, sort_order='bottom'):
        '''
        compactly print the matrix product operators at each site with local operators represented as just symbols
        :param sort_order: the node with the identity edge from the head node is at the "bottom"/"top" of the layer
        :return: None
        '''
        symbol_mpo = self.GenSymbolMPO(sort_order)
        print('')
        for layer_idx in range(0, len(symbol_mpo)):
            # print row-by-row
            for i in range(0, len(symbol_mpo[layer_idx])):
                print(symbol_mpo[layer_idx][i])
            print('')

    def PrintBondDimension(self):
        '''
        print the MPO bond dimensions
        :return: None
        '''
        print('')
        for layer_idx in range(0, self.site_num + 1):
            print(len(self.locator[layer_idx]), end=' ')
        print('')


    def SortLayer(self, sort_order='bottom'):
        '''
        sort each layer to obtain translational symmetry
        :param sort_order: the node with the identity edge from the head node is at the "bottom"/"top" of the layer
        :return: None
        '''
        # define the identity
        opr_shape = self.head.list_edge_out[0].named_opr.data.shape
        Id = named_data(self.name_Id, np.identity(opr_shape[0]))

        if (sort_order == 'top'):
            # 'top' example-TFI chain:
            # [[Id, Sz, Sx],
            #  [0,  0,  Sz],
            #  [0,  0,  Id]]
            last_node = self.head
            while True:
                # search for edge with identity operator
                the_edge = last_node.SearchEdgeOut(Id, 1)
                if (the_edge == None):
                    break
                else:
                    if (the_edge.node_to.bond_idx != 0):
                        # keep the identity edge from the head node at the top
                        self.LayerTop(the_edge.node_to)
                    last_node = the_edge.node_to

            last_node = self.tail
            while True:
                # search for edge with identity operator
                the_edge = last_node.SearchEdgeIn(Id, 1)
                if (the_edge == None):
                    break
                else:
                    if (the_edge.node_from.bond_idx != len(self.locator[the_edge.node_from.layer_idx]) - 1):
                        # keep the identity edge from the head node at the bottom
                        self.LayerBottom(the_edge.node_from)
                    last_node = the_edge.node_from

        elif (sort_order == 'bottom'):
            # 'bottom' example-TFI chain:
            # [[Id,  0,   0],
            #  [Sz,  0,   0],
            #  [Sx,  Sz,  Id]]
            last_node = self.head
            while True:
                # search for edge with identity operator
                the_edge = last_node.SearchEdgeOut(Id, 1)
                if (the_edge == None):
                    break
                else:
                    if (the_edge.node_to.bond_idx != len(self.locator[the_edge.node_to.layer_idx]) - 1):
                        # keep the identity edge from the head node at the bottom
                        self.LayerBottom(the_edge.node_to)
                    last_node = the_edge.node_to

            last_node = self.tail
            while True:
                # search for edge with identity operator
                the_edge = last_node.SearchEdgeIn(Id, 1)
                if (the_edge == None):
                    break
                else:
                    if (the_edge.node_from.bond_idx != 0):
                        # keep the identity edge from the head node at the top
                        self.LayerTop(the_edge.node_from)
                    last_node = the_edge.node_from
        else:
            print('---Error: illegal sort_order!----')
            exit(1)

    def LayerTop(self, node_top):
        '''
        shift the given node to the top of the layer
        :param node_top: the node to be shifted at the top
        :return: None
        '''
        # change bond_idx smaller than that of node_top
        for the_node in self.locator[node_top.layer_idx][0:node_top.bond_idx]:
            the_node.bond_idx += 1
        # delete the record of node_top in the locator
        del self.locator[node_top.layer_idx][node_top.bond_idx]
        # set bond_idx of node_top to be the smallest
        node_top.bond_idx = 0
        # insert node_top to the top of the layer
        self.locator[node_top.layer_idx].insert(0, node_top)

    def LayerBottom(self, node_bot):
        '''
        shift the given node to the bottom of the layer
        :param node_bot: the node to be shifted at the bottom
        :return: None
        '''
        # change bond_idx larger than that of node_bot
        for the_node in self.locator[node_bot.layer_idx][(node_bot.bond_idx + 1):]:
            the_node.bond_idx -= 1
        # delete the record of node_bot in the locator
        del self.locator[node_bot.layer_idx][node_bot.bond_idx]
        # set bond_idx of node_bot to be the largest
        node_bot.bond_idx = len(self.locator[node_bot.layer_idx])
        # insert node_top to the bottom of the layer
        self.locator[node_bot.layer_idx].append(node_bot)