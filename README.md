# Topology-in-Condensed-Matter_code
Public

IIT Bombay, Course Project
<br>
Author- RASHMI MEENA

Code_1 : Band Diagram

    import numpy as np
    import scipy.linalg as la
    import matplotlib.pyplot as plt
    
    # Define the range of k values
    k = np.linspace(-np.pi, np.pi, 100)
    
    # Initialize empty lists to store eigenvalues and k values
    eigenvalues_1 = []
    eigenvalues_2 = []
    k_values = []
    
    # Define the parameters
    t_B = 1
    t_A = 0
    
    # Loop over each k value
    for i in range(len(k)):
    # Create the Hamiltonian matrix A for each k value
    A = np.array([[0, t_B * np.exp(1j * k[i]) + t_A],
                  [t_A + t_B * np.exp(-1j * k[i]), 0]])
    
    # Calculate the eigenvalues
    eig = la.eigvalsh(A)
    
    # Append the eigenvalues to the lists
    eigenvalues_1.append(eig[0].real)
    eigenvalues_2.append(eig[1].real)
    
    # Append the k value to the list
    k_values.append(k[i])

    # Plot the results
    plt.plot(k_values, eigenvalues_1, label='Eigenvalue 1')
    plt.plot(k_values, eigenvalues_2, label='Eigenvalue 2')
    
    plt.xlabel("KX")
    plt.ylabel("Energy")
    plt.title("E vs K plot [-$\pi$,$\pi$], t_A=0, t_B=1", color="darkorange")
    plt.legend()
    plt.grid(True)
    plt.show()
    
Code_2 : Winding the origin for different parameter value

    import matplotlib.pyplot as plt
    import numpy as np
    
    t_A = 0.5
    t_B = 1.0
    k = np.linspace(-np.pi, np.pi, 100)
    X = t_A + t_B * np.cos(k)
    Y = t_B * np.sin(k)
    
    plt.axvline(0, c="black")
    plt.axhline(0, c="black")
    
    plt.plot(X, Y)
    plt.grid(True)
    plt.title("Winding in h plane (t_A=0.5, t_B=1.0)")
    plt.xlabel("h_x")
    plt.ylabel("h_y")
    plt.show()
