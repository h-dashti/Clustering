#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class HK2d :
    '''
    Implementation of Union-Find Algorithm.

    labels array: has the meaning that labels[x] is an
    alias for the label x; by following this chain until
    x == labels[x], you can find the canonical name of an
    equivalence class.	The labels start at one; labels[0] is
    a special value indicating the highest label already used.

    site array: the  label of clusters
    '''


    def __init__(self, shape):
        '''
        	shape: the shape of the lattice
        '''
        self.__size = np.prod(shape)
        self.__shape = shape
        self.__labels = np.empty(self.__size // 2 + 1, dtype=int)
        self.__site = np.empty(self.__size, dtype=int)
        self.__n_clusters = None

    #==============================================

    def __find(self, x):
        '''
        returns the canonical label for the equivalence class containing x
        '''
        y = x
        while self.__labels[y] != y:
            y = self.__labels[y]
        while self.__labels[x] != x:
            z = self.__labels[x]
            self.__labels[x] = y
            x = z

        return y
    #==============================================

    def __union(self, x, y) :
        '''
        joins two equivalence classes and
        returns the canonical label of the resulting class.
        '''
        value = self.__find(y)
        self.__labels[self.__find(x)] = value
        return value
    #==============================================

    def __make_set(self):
        '''
        creates a new equivalence class and returns its label
        '''
        self.__labels[0] += 1
        self.__labels[self.__labels[0]] = self.__labels[0]
        return self.__labels[0]
    #==============================================

    def __impose_pbc(self, axes):
        '''
        To consider periodic boundary conditon along axes
        '''
        if axes[0]:
            #  x1=0,  x2=lx-1
            for i1 in range(self.__shape[1]):
                first = self.__site[i1]
                if first == 0:  continue
                second = self.__site[i1 + self.__size - self.__shape[1]]
                if second == 0: continue
                self.__union(first, second)

        if axes[1]:
            #  y1=0,  y2=ly-1
            for i1 in range(0, self.__size, self.__shape[1]):
                first = self.__site[i1]
                if first == 0:  continue
                second = self.__site[i1 + self.__shape[1] - 1]
                if second == 0: continue
                self.__union(first, second)

    #==============================================

    def clustering(self, occupied,
                   trgt = 1,
                   pbc_axes = [False, False]
                   ):
        '''
        Labeling the clusters in occupied array.
        The pbc_axes is periodic boundary conditons along x and y axes.
        '''

        # first set the site array
        if occupied.ndim == 1:
            self.__site = np.where(occupied == trgt, 1, 0)
        else:
            self.__site = np.where(occupied.ravel() == trgt, 1, 0)

        self.__labels[0] = 0

        size = np.prod(self.__shape)
        l1 = self.__shape[1]

        for i in range(size):
            if self.__site[i] == 0:
                continue

            x = i // l1         # i = x*l1+y
            y = i - x * l1

            up = self.__site[i - l1]  if x > 0  else 0 #  (x-1, y)
            left = self.__site[i - 1] if y > 0  else 0 #  (x, y-1)

            if up == 0 and left == 0:  # neither a label above nor left
                self.__site[i] = self.__make_set() # make a new cluster

            elif up == 0 and left != 0:   # one neighbour down
                self.__site[i] = left
            elif up != 0 and left == 0:   # one neighbour left
                self.__site[i] = up
            else:
                self.__site[i] = self.__union(up, left)
        # end for i

        ##############
        # addiing PBC
        ##############
        self.__impose_pbc(pbc_axes)


        # We apply the relabeling to the site array.
        # we create a mapping from the canonical labels
		# determined by union/find into a new set of canonical labels,
        # which are guaranteed to be sequential.
        new_labels = np.zeros(self.__labels.shape, dtype=int)

        for i in range(size):
            if self.__site[i] == 0:
                continue
            x = self.__find(self.__site[i])
            if new_labels[x] == 0:
                new_labels[0] += 1;
                new_labels[x] = new_labels[0]
            self.__site[i] = new_labels[x]

        self.__n_clusters = new_labels[0]
        #return self.__site.reshape(occupied.shape), new_labels[0]

    #==============================================
    @property
    def site(self):
        '''
        Return site array wich is a row wise one
        '''
        return self.__site

    @property
    def n_clusters(self):
        return self.__n_clusters

    @property
    def shape(self):
        return self.__shape

#################################################################3


if __name__ == '__main__':

    l = 5
    #np.random.seed(3345)
    arr = [ 1 if np.random.random() < 0.5 else 0 for i in range(l**2)]
    arr = np.int8(arr).reshape(l, l)

    hk = HK2d(arr.shape)

    hk.clustering(arr, pbc_axes=[True,True])
    print(hk.site.reshape(arr.shape))
    print(hk.n_clusters)







