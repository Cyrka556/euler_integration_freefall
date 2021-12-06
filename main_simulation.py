#Definitions
import math
import matplotlib.pyplot as plt
import numpy as np



#Constant
g = 9.81 #m/s^2


#Function definitions
def spharea(r):
    """User only requires an input of the sphere's radius to use surface area"""

    area = math.pi*r**2
    return area




def drag(Cd,A,p):
    """Assumes a sphere, constant Cd(0.47) and area based on inputted radius"""
    k = (Cd*A*p)/2
    return k




def height(yo, m, k, t):
    """Height of ball after being dropped from initial height (yo), based eqn 8"""
    height = (yo - ((m/k)*np.log(np.cosh((np.sqrt((k*g)/m))*t))))
    return height




def tmax(k,m,yo):
    maximum = (math.sqrt(m/(k*g)))*(math.acosh(math.exp((k/m)*yo)))

    return maximum




def speed(m,k,t):
    """Speed from equation 9, not negative as Speed was what's calculated, not the direction"""
    coef = np.sqrt((m*g)/k)
    inside = np.sqrt(((k*g)/m)*t)
    tanh = np.tanh(inside)
    velocity = (-1)*coef*tanh

    return velocity







def euler(m,p,A,Cd,dt,yo):
    """Euler integration method, stops at y = 0, with varying density included"""

    #Assumption that the area remains constant, doesn't change orientation (not very accurate)


    #Calculted values
    k = drag(Cd, A, p)

    time = [0]
    speed = [0]
    dist = [yo]

    vy = 0
    y = yo
    t=0

    while y >= 0:

        t += dt
        time.append(t)

        vy = (speed[-1] - dt * (g + ((k / m) * (np.abs(speed[-1])) * (speed[-1]))))
        speed.append(vy)

        y = (dist[-1] + dt*(speed[-1]))
        dist.append(y)



    return time, speed, dist





def eulerp():
    """Euler integration method, stops at y = 0, with varying density included"""

    #List of constants for Felix jump from 37640m, time step of 0.1s
    #Assumption that the area remains constant, doesn't change orientation (not very accurate)
    m = 110 #kg
    po = 1.2 #kg/m^3
    A = 2 #m^2
    Cd = 1
    dt = 0.1
    yo = 37640
    h = 7.64 #km^2


    time = [0]
    speed = [0]
    dist = [yo]
    accelaration = [g]


    vy = 0
    y = yo
    t=0

    while y >= 0:
        yk = y/1000
        p = po*(math.exp((-yk)/h))

        k = drag(Cd,A,p)

        t += dt
        time.append(t)

        vy = (speed[-1] - dt * (g + ((k / m) * (np.abs(speed[-1])) * (speed[-1]))))
        speed.append(vy)

        y = (dist[-1] + dt*(speed[-1]))
        dist.append(y)

        acc = -(g + ((k / m) * (np.abs(speed[-1])) * (speed[-1])))
        accelaration.append(acc)


    return time, speed, dist, accelaration





def ModifiedEuler():
    """Uses 2 gradients around a point to more precisely measure the gradient at the midpoint, becoming more accurate than the normal euler method"""


    dt = 0.1  # s



    time, speed, dist, acc = eulerp()

    ModSpeed = [0]
    ModDist = [37640]



    for i in range(0,len(speed)-1):

        yn = ModDist[i] + (dt/2)*(speed[i] + speed[i+1])

        ModDist.append(yn)

        vn = ModSpeed[i] + (dt/2)*(acc[i] + acc[i+1])
        ModSpeed.append(vn)




    HeightDiff = []
    SpeedDiff = []
    for i in range(0, len(ModDist)):
        heightdiff = np.abs(ModDist[i] - dist[i])
        HeightDiff.append(heightdiff)

        speeddiff = -np.abs(speed[i] - ModSpeed[i])
        SpeedDiff.append(speeddiff)




    return ModDist, HeightDiff, ModSpeed, SpeedDiff







#Main menu to navigate different methods

loop = 5
while loop<6:
    section = input('Go to\na - Equation Method\nb - Euler Integration Method\nc - Euler Integration with varying density\nd - Modified Euler Integration with Varying Density\ne - Comparison of Methods\nq - Quit\n')
    print()


    if section == 'a':
        aloop = 5
        while aloop<6:
            print('This simulates a ball falling with drag coefficient 0.47 and constant air density (1.2 kg/m^3)')
            print()



            m = float(input('Mass of the sphere (kg) = '))
            radius = float(input('Radius of the sphere (m) = '))
            yo = float(input('Initial height (m) = '))

            A = spharea(radius)

            k = drag(0.47, A, 1.2)

            maxt = tmax(k, m, yo)

            time = np.linspace(0, maxt, (maxt * 1000))
            y = height(yo, m, k, time)

            points = int(maxt * 1000)
            time = np.linspace(0, maxt, points)
            vy = np.zeros(points)
            for i in range(points):
                vy[i] = speed(m, k, i)

            plt.figure(1)
            plt.subplot(311)
            plt.plot(time, y, 'g')
            plt.ylabel('Displacement (m)')
            plt.xlabel('Time (s)')

            plt.subplot(312)
            plt.subplots_adjust(hspace=0.8)
            plt.plot(time, vy, 'b')
            plt.ylabel('Speed (m/s)')
            plt.xlabel('Time (s)')

            plt.show()


            print()
            option = input('Choose an option\nc - Rerun this option\nr - Return to main menu\n')
            print()

            if option == 'c':
                pass

            elif option == 'r':
                break




    elif section == 'b':
        bloop = 5
        while bloop<6:
            print('Simulation of the jump from 37640m\nMass = 110kg\nDrag coefficient = 1\nConstant density of 1.2 kg/m^3\nArea = 2 m^2\nIntegration step of 0.1 s')
            print()

            plt.figure(1)
            plt.subplot(211)
            plt.plot((euler(110, 1.2, 2, 1, 0.1, 37640))[0], (euler(110, 1.2, 2, 1, 0.1, 37640))[2], 'r')
            plt.ylabel('Distance (m)')
            plt.xlabel('Time (s)')

            plt.subplot(212)
            plt.plot((euler(110, 1.2, 2, 1, 0.1, 37640))[0], (euler(110, 1.2, 2, 1, 0.1, 37640))[1])
            plt.ylabel('Velocity (m/s)')
            plt.xlabel('Time (s)')

            plt.show()

            print()
            option = input('Choose an option\nc - Rerun this option\nr - Return to main menu\n')
            print()


            if option == 'c':
                pass

            elif option == 'r':
                break



    elif section == 'c':
        print('Simulation of the jump from 37640m\nMass = 110kg\nDrag coefficient = 1\nVarying density\nArea = 2 m^2\nIntegration step of 0.1 s')
        print()


        cloop = 5
        while cloop<6:
            plt.figure(1)
            plt.subplot(211)
            plt.subplots_adjust(hspace=0.29)
            plt.plot((eulerp())[0], (eulerp())[2], 'r')
            plt.ylabel('Distance (m)')
            plt.xlabel('Time (s)')

            plt.subplot(212)
            plt.plot((eulerp())[0], (eulerp())[1])
            plt.ylabel('Velocity (m/s)')
            plt.xlabel('Time (s)')

            plt.show()

            print()
            option = input('Choose an option\nc - Rerun this option\nr - Return to main menu\n')
            print()

            if option == 'c':
                pass

            elif option == 'r':
                break

    elif section == 'd':
        print('Simulation of the jump from 37640m\nMass = 110kg\nDrag coefficient = 1\nVarying density\nArea = 2 m^2\nIntegration step of 0.1 s, under modified euler method')
        print()


        dloop = 5
        while dloop < 6:
            fig1 = plt.figure(1)
            fig1.subplots_adjust(hspace=0.35)
            plt.subplot(211)
            plt.plot(eulerp()[0], ModifiedEuler()[0], 'r')
            plt.ylabel('Distance (m)')
            plt.xlabel('Time (s)')


            plt.subplot(212)
            plt.plot(eulerp()[0], ModifiedEuler()[2], 'b')
            plt.ylabel('Velocity (m/s)')
            plt.xlabel('Time (s)')

            plt.show()

            print()
            option = input('Choose an option\nc - Rerun this option\nr - Return to main menu\n')
            print()

            if option == 'c':
                pass

            elif option == 'r':
                break



    elif section == 'e':

        print('Comparison of methods, choose from list below')
        print()

        comparison = input('Choose from:\na - Functions and euler step for constant density\nb - Euler and Modified Euler for varying density\nr - Return to the main menu\n')

        if comparison == 'a':
            A = spharea(0.2)

            k = drag(0.47, A, 1.2)

            tmaximum = tmax(k, 0.5, 1000)

            acloop = 5
            while acloop <6:
                print('Comparison for function and euler methods for travelling of a sphere with:\nMass = 0.5kg\nCd = 0.47\nConstant Density = 1.2 kg/m^3\nyo = 1000m\n')


                time = np.linspace(0, tmaximum, (tmaximum * 1000))
                y = height(1000, 0.5, k, time)


                points = int(tmaximum * 1000)
                time = np.linspace(0, tmaximum, points)
                vy = np.zeros(points)
                for i in range(points):
                    vy[i] = speed(0.5, k, i)


                plt.figure(1)
                plt.subplot(221)
                plt.subplots_adjust(bottom=0.02, top=0.78, wspace=0.4, hspace=0.01)
                plt.plot((euler(0.5, 1.2, spharea(0.2), 0.47, 0.1, 1000))[0], (euler(0.5, 1.2, spharea(0.2), 0.47, 0.1, 1000))[2], 'r', time, y, 'b')
                plt.ylabel('Distance (m)')
                plt.xlabel('Time (s)')
                plt.legend(['Euler Method', 'Function'], loc='center right',bbox_to_anchor=(0.975, 1.07))



                plt.subplot(222)
                plt.plot((euler(0.5, 1.2, spharea(0.2), 0.47, 0.1, 1000))[0], (euler(0.5, 1.2, spharea(0.2), 0.47, 0.1, 1000))[1], 'r', time, vy, 'b')
                plt. ylabel('Velocity (m/s)')
                plt.xlabel('Time (s)')


                plt.show()


                print()
                option = input('Choose an option\nc - Rerun this option\nr - Return to main menu\n')
                print()

                if option == 'c':
                    pass

                elif option == 'r':
                    break


        elif comparison == 'b':
            bcloop = 5
            while bcloop < 6:

                print('Comparison of euler and modified euler method for Baumgartner jump\nMass = 110kg\nVarying Density with po = 1.2 kg/m^3\n'
                      'Area = 2 m^2\nyo = 37640\n')


                fig1 = plt.figure(1)
                fig1.subplots_adjust(top=0.93, hspace=0.35)
                plt.subplot(211)
                plt.plot(eulerp()[0], ModifiedEuler()[0], 'r', eulerp()[0], eulerp()[2], 'purple', eulerp()[0],
                         ModifiedEuler()[1], 'darkblue', linewidth = 0.7)
                plt.ylabel('Distance (m)')
                plt.xlabel('Time (s)')
                plt.legend(['Modified Euler', 'Euler', 'Difference'], loc='center',
                           bbox_to_anchor=(0.975, 0.5))

                plt.subplot(212)
                plt.plot(eulerp()[0], ModifiedEuler()[2], 'r', eulerp()[0], eulerp()[1], 'purple', eulerp()[0],
                         ModifiedEuler()[3], 'darkblue', linewidth = 0.7)
                plt.ylabel('Velocity (m/s)')
                plt.xlabel('Time (s)')

                plt.show()

                print()
                option = input('Choose an option\nc - Rerun this option\nr - Return to main menu\n')
                print()

                if option == 'c':
                    pass

                elif option == 'r':
                    break


        elif comparison == 'r':
            print('I will return you to the main menu')
            print()
            break

        else:
            print('I will return you to the main menu')
            print()
            break


    elif section == 'q':
        print('Ending program')
        break


    else:
        print('I will assume you want to end the program')
        break
