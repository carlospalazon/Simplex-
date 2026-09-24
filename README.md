# Simplex-

## Descripció del projecte
Aquest repositori conté la implementació completa i robusta de l'algorisme Símplex per a la resolució de problemes de Programació Lineal (PL). Desenvolupat íntegrament en Python a través de la classe `SimplexSolver`, l'algorisme està dissenyat per resoldre de principi a fi casos d'optimització complexos, detectant i interpretant escenaris amb solucions òptimes, problemes infactibles i problemes no acotats.

Per validar-ne la fiabilitat, el model s'ha sotmès a proves sobre dos conjunts de dades algorítmics (Conjunt 34 i Conjunt 37), resolent amb èxit un total de 8 problemes d'optimització amb diferents topologies.

## Característiques principals
- **Mètode de Dues Fases:** Garanteix l'arrencada de l'algorisme en qualsevol escenari, independentment de les restriccions inicials, introduint i gestionant variables artificials.
- **Regla de Bland:** Implementació estricta en el pivoteig per prevenir el ciclatge infinit en problemes degenerats, assegurant la convergència finita del mètode.
- **Actualització Eficient de la Base:** Utilització de matrius Eta per actualitzar la inversa de la base $B^{-1}$, la qual cosa redueix el cost computacional a $O(m^2)$ per iteració en lloc de calcular la inversa des de zero.
- **Gestió Activa de la Degeneració:** Tractament segur de variables artificials que romanen a la base amb valor zero en finalitzar la Fase I, mantenint-les de forma innòcua a la Fase II.
