**********************
**Chemical reactions**
**********************

This extension enables the rate-controlled stochastic curing of polymer
systems, either for chain growth of step growth, depending on the
parameters.

The variables typeA, typeB, minStateA, minStateB, maxStateA, maxStateB
control the particles that enter the curing reaction

.. math::
  A^a + B^b \\rightarrow A^{a+deltaA}-B^{b+deltaB}

where A and B may possess additional bonds not shown.

An extra bond is added between A and B whenever the state of A and B falls
into the defined range by variables min/max state.
The conditions are as follow:

.. math::

   a >= minStateA \\land stateA < maxStateA

the same holds for the particle B. Both conditions for *A* and *B* have to match.
The second condition is probabilistic:

.. math::

   k\Delta t \Phi < W

where :math:`k` is a *kinetic rate*, :math:`\Delta t` is an integrator
timestep, :math:`\Phi` is an interval and :math:`W` is a random number from
uniform distribution.

In addition if the intramolecular property is set to true (by default) then
the reaction only happen between heterogeneous molecules.

The reaction proceeds by testing for all possible (A,B) pairs and
selects them only at a given rate. It works in parallel, by gathering
first the successful pairs between neighboring CPUs and ensuring that
each particle enters only in one new bond per reaction step.

Example
######################

**Creating the integrator extension**

>>> ar = espressopp.integrator.ChemicalReaction(
>>>     system,
>>>     verletList,
>>>     system.storage,
>>>     interval)

**Creates synthesis reaction**

>>> r_type_1 = espressopp.integrator.Reaction(
>>>     type_a=ar_type_M,
>>>     type_b=ar_type_B,
>>>     delta_a=1,
>>>     delta_b=1,
>>>     min_state_a=0,
>>>     min_state_b=0,
>>>     max_state_a=2,
>>>     max_state_b=1,
>>>     rate=1000,
>>>     cutoff=ar_cutoff,
>>>     fpl=fpl_a_b
>>>     )

Add the reaction to the integrator extension

>>> ar.add_reaction(r_type_1)

Add the extension to the integrator

>>> integrator.addExtension(ar)


.. automodule:: espressopp.integrator.ChemicalReaction
   :members:
