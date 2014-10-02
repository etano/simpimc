import IlkkaSquarer
import DavidSquarer

def run(units,particles,potential,squarer,breakup,objects):

    if squarer['type'] is 'Ilkka':

        print '\n**** Performing squaring ****\n'
        print 'WARNING: Assuming 3D Coulomb interaction'
        IlkkaSquarer.Square(particles,squarer)

        print '\n**** Performing breakup ****\n'
        IlkkaSquarer.Breakup(particles,potential,squarer,breakup,objects)

    elif squarer['type'] is 'David':

        print '\n**** Performing breakup ****\n'
        DavidSquarer.Breakup(units,particles,potential,squarer,breakup,objects)

        print '\n**** Performing squaring ****\n'
        DavidSquarer.Square(particles,squarer,objects)
