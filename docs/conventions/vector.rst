.. _conventions_vector:

Vector representation
=====================

In GooseFEM there are three ways to represent vectors. In particular, a vector field (e.g. the displacement) can be collected:

* per node (denoted "nodevec", :ref:`see below <conventions_storage>`),
* per degree-of-freedom (denoted "dofval", :ref:`see below <conventions_storage>`),
* per element  (denoted "elemvec", :ref:`see below <conventions_storage>`).

.. warning::

  Watch out with the conversion from one representation to the other as downsizing can be done in more than one way, see :ref:`conventions_vector_conversion`.

Consider a simple two-dimensional mesh of just two elements, and a displacement vector per node:

.. image:: figures/data-representation.svg
  :width: 400px
  :align: center

|

Collected per node (nodevec)
----------------------------

.. math::

  \texttt{disp} =
  \begin{bmatrix}
    u_x^{(0)} & u_y^{(0)} \\
    u_x^{(1)} & u_y^{(1)} \\
    u_x^{(2)} & u_y^{(2)} \\
    u_x^{(3)} & u_y^{(3)} \\
    u_x^{(4)} & u_y^{(4)} \\
    u_x^{(5)} & u_y^{(5)}
  \end{bmatrix}

Collected per degree-of-freedom (dofval)
----------------------------------------

The following definition

.. math::

  \texttt{dofs} =
  \begin{bmatrix}
     0 &  1 \\
     2 &  3 \\
     4 &  5 \\
     6 &  7 \\
     8 &  9 \\
    10 & 11
  \end{bmatrix}

gives:

.. math::

  \texttt{u} =
  \big[
    u_x^{(0)} \,
    u_y^{(0)} \,
    u_x^{(1)} \,
    u_y^{(1)} \,
    u_x^{(2)} \,
    u_y^{(2)} \,
    u_x^{(3)} \,
    u_y^{(3)} \,
    u_x^{(4)} \,
    u_y^{(4)} \,
    u_x^{(5)} \,
    u_y^{(5)}
  \big]^T

Whereby "dofs" can be used to:

* **Reorder** "u" such that is can be easily (even directly) partitioned. For example, consider that all :math:`x`-coordinates are *Prescribed* and all :math:`y`-coordinates are *Unknown*. In particular,

  .. math::

    \texttt{dofs} =
    \begin{bmatrix}
       6 & 0 \\
       7 & 1 \\
       8 & 2 \\
       9 & 3 \\
      10 & 4 \\
      11 & 5
    \end{bmatrix}

  gives

  .. math::

    \texttt{u} =
    \big[
      u_y^{(0)} \,
      u_y^{(1)} \,
      u_y^{(2)} \,
      u_y^{(3)} \,
      u_y^{(4)} \,
      u_y^{(5)} \, \;
      u_x^{(0)} \,
      u_x^{(1)} \,
      u_x^{(2)} \,
      u_x^{(3)} \,
      u_x^{(4)} \,
      u_x^{(5)}
    \big]^T
    =
    \big[
      \texttt{u}_u \, \;
      \texttt{u}_p
    \big]^T

  which allows

  .. math::

    \texttt{u}_u &= \texttt{u[:6]} \\
    \texttt{u}_p &= \texttt{u[6:]}

  |

* **Eliminate** dependent nodes. For example, suppose that the displacement of all top nodes is equal to that of the bottom nodes. In this one could:

  .. math::

    \texttt{dofs} =
    \begin{bmatrix}
       0 & 1 \\
       2 & 3 \\
       4 & 5 \\
       0 & 1 \\
       2 & 3 \\
       4 & 5
    \end{bmatrix}
    \qquad
    \rightarrow
    \qquad
    \texttt{u} =
    \begin{bmatrix}
      u_0 \\
      u_1 \\
      u_2 \\
      u_3 \\
      u_4 \\
      u_5
    \end{bmatrix}
    \quad
    \leftrightarrow
    \quad
    \texttt{disp} =
    \begin{bmatrix}
      u_0 & u_1 \\
      u_2 & u_3 \\
      u_4 & u_5 \\
      u_0 & u_1 \\
      u_2 & u_3 \\
      u_4 & u_5
    \end{bmatrix}

.. note::

  :ref:`Vector` applies the reordering itself. One does not need to change "dofs", but one simply supplies "iip".

Collected per element (elemvec)
-------------------------------

For this example:

.. math::

  \texttt{conn} =
  \begin{bmatrix}
    0 & 1 & 4 & 3 \\
    1 & 2 & 5 & 4
  \end{bmatrix}

The storage per node proceeds in

.. math::

  \texttt{shape(ue)}
  &= \left[ n_\text{elements} \times n_\text{nodes-per-element} \times n_\text{dim} \right]
  \\
  &= \left[ 2 \times 4 \times 2 \right]

In particular:

.. math::

  \texttt{ue[0,:,:]} =
  \begin{bmatrix}
    u_x^{(0)} & u_y^{(0)} \\
    u_x^{(1)} & u_y^{(1)} \\
    u_x^{(4)} & u_y^{(4)} \\
    u_x^{(3)} & u_y^{(3)} \\
  \end{bmatrix}

and

.. math::

  \texttt{ue[1,:,:]} =
  \begin{bmatrix}
    u_x^{(1)} & u_y^{(1)} \\
    u_x^{(2)} & u_y^{(2)} \\
    u_x^{(5)} & u_y^{(5)} \\
    u_x^{(4)} & u_y^{(4)} \\
  \end{bmatrix}

.. _conventions_vector_conversion:

Conversion
----------

Conversion to a larger representation (upsizing) can always be done uniquely, however, conversion to a more compact representation (downsizing) can be done in two ways. In particular:

+---------+---------+-------------------+------------------------------+
| From    | To      | Function          | Remarks                      |
+=========+=========+===================+==============================+
| dofval  | nodevec | asNode(...)       | unique                       |
+---------+---------+-------------------+------------------------------+
| dofval  | elemvec | asElement(...)    | unique                       |
+---------+---------+-------------------+------------------------------+
| nodevec | elemvec | asElement(...)    | unique                       |
+---------+---------+-------------------+------------------------------+
| nodevec | dofval  | asDofs(...)       | overwrites reoccurring items |
+---------+---------+-------------------+------------------------------+
| elemvec | dofval  | asDofs(...)       | overwrites reoccurring items |
+---------+---------+-------------------+------------------------------+
| elemvec | nodevec | asNode(...)       | overwrites reoccurring items |
+---------+---------+-------------------+------------------------------+
| nodevec | dofval  | assembleDofs(...) | adds reoccurring items       |
+---------+---------+-------------------+------------------------------+
| elemvec | dofval  | assembleDofs(...) | adds reoccurring items       |
+---------+---------+-------------------+------------------------------+
| elemvec | nodevec | assembleNode(...) | adds reoccurring items       |
+---------+---------+-------------------+------------------------------+
