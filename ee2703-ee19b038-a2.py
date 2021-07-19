"""
        EE2703:Applied Programming Lab
        Assignment 2: Spice - Part 2
        NAME: Manikandan Sritharan
        ROLL NO: EE19B038
"""

"""
        INPUT FORMAT: python3 <name_of_python_file> "<path_of_the_netlist_file>"
        (Example: python3 ee2703_ee19b038_a2.py "ckt1.netlist")
        Note: path_of_the_netlist_file can be the full path or relative path.
"""

# Necessary Libraries
import sys
import numpy as np
import pandas as pd
import cmath

# For easier control of various symbols, we define them here.
CIRCUIT = '.circuit'
END = '.end'
R = "R"
C = "C"
I = "L"
IVS = "V"
ICS = "I"
VCVS = "E"
VCCS = "G"
CCVS = "H"
CCCS = "F"

# Defining classes to model various Circuit Elements
# ni's represents nodes
class resistor:
    def __init__(self, name, n1, n2, value):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.value = float(value)
        

class inductor:
    def __init__(self, name, n1, n2, value):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.value = float(value)
        

class capacitor:
    def __init__(self, name, n1, n2, value):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.value = float(value)
        

class voltageSource:
    def __init__(self, name, n1, n2, value, phase=0):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.value = float(value)
        self.phase = float(phase)

class currentSource:
    def __init__(self, name, n1, n2, value, phase=0):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.value = float(value)
        self.phase = float(phase)

class vcvs:
    def __init__(self, name, n1, n2, n3, n4, value):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        self.n4 = n4
        self.value = float(value)
        

class vccs:
    def __init__(self, name, n1, n2, n3, n4, value):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        self.n4 = n4
        self.value = float(value)
        

class ccvs:
    def __init__(self, name, n1, n2, Vsource, value):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.Vsource = Vsource
        self.value = float(value)

class cccs:
    def __init__(self, name, n1, n2, Vsource, value):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.Vsource = Vsource
        self.value = float(value)

argumentlist = sys.argv

# Checking if the command line statement is correct
if(len(argumentlist) != 2):
    print("Invalid number of command line arguments\n")
    sys.exit()

filename = argumentlist[1]

if(not filename.endswith(".netlist")):
    print("Not a Netlist.\n")
    sys.exit()
try:
    f = open(filename)
except Exception:
    print("Enter valid file address\n")
    sys.exit()


# Reading the data from the file
datalines = f.readlines()
netlist_data = []

#Removing comments in the circuit definition
PI = np.pi
for line in datalines:
    netlist_data.append(line.split('#')[0].split('\n')[0])
    # if .ac line is present, we need to get the frequency
    freq = 1e-50 
    if(line[:3] == '.ac'):
        freq = float(line.split()[2])
        w = 2*PI*freq

# Determining the section of the code that contains the circuit definition. (if it exists)
try:
    start_index = netlist_data.index(CIRCUIT)
except ValueError:
    print("The line "+CIRCUIT+" is not present in the netlist. Hence invalid netlist.\n")
    sys.exit()

try:
    end_index = netlist_data.index(END)
except ValueError:
    print("The line "+END+" is not present in the netlist. Hence invalid netlist.\n")
    sys.exit()

circuit_data = netlist_data[start_index+1:end_index]

#Closing file as required data has been extracted
f.close()

# List of Nodes and set of Circuit Elements
nodes = []
cirElements = {R: [], C: [], I: [], IVS: [], ICS: [], VCVS: [], VCCS: [], CCVS: [], CCCS: []}

for line in circuit_data:
    tokens = line.split()
    #print(*tokens)
    # Adding new nodes to the list (if any)
    if tokens[1] not in nodes:
        nodes.append(tokens[1])
    elif tokens[2] not in nodes:
        nodes.append(tokens[2])
    
    #Adding the circuit element to cirElements
    # Resistor
    if tokens[0][0] == R:
        cirElements[R].append(resistor(tokens[0], tokens[1], tokens[2], tokens[3]))
        
    # Capacitor
    elif tokens[0][0] == C:
        cirElements[C].append(capacitor(tokens[0], tokens[1], tokens[2], tokens[3]))
        
    # Inductor
    elif tokens[0][0] == I:
        cirElements[I].append(inductor(tokens[0], tokens[1], tokens[2], tokens[3]))
        
    # Voltage Source
    elif tokens[0][0] == IVS:
        if len(tokens) == 5: # DC Source
            cirElements[IVS].append(voltageSource(tokens[0], tokens[1], tokens[2], float(tokens[4])))
        elif len(tokens) == 6: # AC Source
            if freq == 1e-50:
                print("Frequency of AC Source is not mentioned.\n")
                sys.exit()
            cirElements[IVS].append(voltageSource(tokens[0], tokens[1], tokens[2], float(tokens[4])/2, tokens[5]))
        
    # Current Source
    elif tokens[0][0] == ICS:
        # DC
        if len(tokens) == 5:
            cirElements[ICS].append(currentSource(tokens[0], tokens[1], tokens[2], float(tokens[4])))
        # AC
        elif len(tokens) == 6:
            if freq == 1e-50:
                print("Frequency of AC Source is not mentioned.\n")
                sys.exit()
            cirElements[ICS].append(currentSource(tokens[0], tokens[1], tokens[2], float(tokens[4])/2, tokens[5]))
        
    # Voltage Controlled Voltage Source
    elif tokens[0][0] == VCVS:
        cirElements[VCVS].append(vcvs(tokens[0], tokens[1], tokens[2], tokens[3], tokens[4], tokens[5]))
        
    # Voltage Controlled Current Source
    elif tokens[0][0] == VCCS:
        cirElements[VCCS].append(vcvs(tokens[0], tokens[1], tokens[2], tokens[3], tokens[4], tokens[5]))
        
    # Current Controlled Voltage Source
    elif tokens[0][0] == CCVS:
        cirElements[CCVS].append(ccvs(tokens[0], tokens[1], tokens[2], tokens[3], tokens[4]))
        
    # Current Controlled Current Source
    elif tokens[0][0] == CCCS:
        cirElements[CCCS].append(cccs(tokens[0], tokens[1], tokens[2], tokens[3], tokens[4]))
        
    else:
        print("Un-identified Circuit Element.\n")
        sys.exit()

#for i in range(len(cirElements[R])):
#    print(cirElements[R][i].name)

# Checking for 'GND', if it exists, placing it at the start of the list of nodes
try:
    nodes.remove('GND')
    nodes = ['GND'] + nodes
except:
    print("No GND node.\n")
    sys.exit()

#print(*nodes)

# Assigning numbers to the nodes and storing it in a dictionary
num_nodes = len(nodes)
numbered_nodes = {nodes[i]:i for i in range(num_nodes)}

#print(*numbered_nodes)

# Number of Voltage Sources
num_VS = len(cirElements[IVS])+len(cirElements[VCVS])+len(cirElements[CCVS])

# Matrices M and b
M = np.zeros((num_nodes+num_VS, num_nodes+num_VS), np.complex)
b = np.zeros((num_nodes+num_VS,), np.complex)

# GND equation
M[0][0] = 1.0

# Resistor equation
for r in cirElements[R]:
    #print(r.n1, r.n2, end='\n')
    if r.n1 != 'GND':
        M[numbered_nodes[r.n1]][numbered_nodes[r.n1]] += 1/r.value
        M[numbered_nodes[r.n1]][numbered_nodes[r.n2]] -= 1/r.value
    if r.n2 != 'GND':
        M[numbered_nodes[r.n2]][numbered_nodes[r.n1]] -= 1/r.value
        M[numbered_nodes[r.n2]][numbered_nodes[r.n2]] += 1/r.value

# Capacitor equation
for c in cirElements[C]:
    if c.n1 != 'GND':
        M[numbered_nodes[c.n1]][numbered_nodes[c.n1]] += complex(0, w*c.value)
        M[numbered_nodes[c.n1]][numbered_nodes[c.n2]] -= complex(0, w*c.value)
    if c.n2 != 'GND':
        M[numbered_nodes[c.n2]][numbered_nodes[c.n1]] -= complex(0, w*c.value)
        M[numbered_nodes[c.n2]][numbered_nodes[c.n2]] += complex(0, w*c.value)

# Inductor equation
for l in cirElements[I]:
    if l.n1 != 'GND':
        M[numbered_nodes[l.n1]][numbered_nodes[l.n1]] += complex(0, -1.0/(w*l.value))
        M[numbered_nodes[l.n1]][numbered_nodes[l.n2]] -= complex(0, -1.0/(w*l.value))
    if l.n2 != 'GND':
        M[numbered_nodes[l.n2]][numbered_nodes[l.n1]] -= complex(0, -1.0/(w*l.value))
        M[numbered_nodes[l.n2]][numbered_nodes[l.n2]] += complex(0, -1.0/(w*l.value))

# Voltage Source equation
for i in range(len(cirElements[IVS])):
    # current through IVS equation
    if cirElements[IVS][i].n1 != 'GND':
        M[numbered_nodes[cirElements[IVS][i].n1]][num_nodes+i] = 1.0
    if cirElements[IVS][i].n2 != 'GND':
        M[numbered_nodes[cirElements[IVS][i].n2]][num_nodes+i] = -1.0
    # Auxiliary equations
    M[num_nodes+i][numbered_nodes[cirElements[IVS][i].n1]] = -1.0
    M[num_nodes+i][numbered_nodes[cirElements[IVS][i].n2]] = +1.0
    b[num_nodes+i] = cmath.rect(cirElements[IVS][i].value, cirElements[IVS][i].phase*PI/180)

# Current Source equation
for i in cirElements[ICS]:
    if i.n1 != 'GND':
        b[numbered_nodes[i.n1]] = -1*i.value
    if i.n2 != 'GND':
        b[numbered_nodes[i.n2]] = i.value

# VCVS equations
for i in range(len(cirElements[VCVS])):
    # current throuh VCVS equation
    if cirElements[VCVS][i].n1 != 'GND':
        M[numbered_nodes[cirElements[VCVS][i].n1]][num_nodes+len(cirElements[IVS])+i] = 1.0
    if cirElements[VCVS][i].n2 != 'GND':
        M[numbered_nodes[cirElements[VCVS][i].n2]][num_nodes+len(cirElements[IVS])+i] = -1.0
    M[num_nodes+len(cirElements[IVS])+i][numbered_nodes[cirElements[VCVS][i].n1]] = 1.0
    M[num_nodes+len(cirElements[IVS])+i][numbered_nodes[cirElements[VCVS][i].n2]] = -1.0
    M[num_nodes+len(cirElements[IVS])+i][numbered_nodes[cirElements[VCVS][i].n3]] = -1.0*cirElements[VCVS][i].value
    M[num_nodes+len(cirElements[IVS])+i][numbered_nodes[cirElements[VCVS][i].n4]] = 1.0*cirElements[VCVS][i].value

# CCVS equations
for i in range(len(cirElements[CCVS])):
    # current through CCVS equation
    if cirElements[VCVS][i].n1 != 'GND':
        M[numbered_nodes[cirElements[CCVS][i].n1]][num_nodes+len(cirElements[IVS])+len(cirElements[VCVS])+i] = 1.0
    if cirElements[VCVS][i].n2 != 'GND':
        M[numbered_nodes[cirElements[VCVS][i].n2]][num_nodes+len(cirElements[IVS])+len(cirElements[VCVS])+i] = -1.0
    M[num_nodes+len(cirElements[IVS])+len(cirElements[VCVS])+i][numbered_nodes[cirElements[CCVS][i].n1]] = 1.0
    M[num_nodes+len(cirElements[IVS])+len(cirElements[VCVS])+i][numbered_nodes[cirElements[CCVS][i].n2]] = -1.0
    M[num_nodes+len(cirElements[IVS])+len(cirElements[VCVS])+i][num_nodes+len(cirElements[IVS])+len(cirElements[VCVS])+i] = -1.0*cirElements[CCVS][i].value

# VCCS equations
for vccs in cirElements[VCCS]:
    if vccs.n1 != 'GND':
        M[numbered_nodes[vccs.n1]][numbered_nodes[vccs.n4]]+=vccs.value
        M[numbered_nodes[vccs.n1]][numbered_nodes[vccs.n3]]-=vccs.value
    if vccs.n2 != 'GND':
        M[numbered_nodes[vccs.n2]][numbered_nodes[vccs.n4]]-=vccs.value
        M[numbered_nodes[vccs.n3]][numbered_nodes[vccs.n3]]+=vccs.value

# CCCS equations
for cccs in cirElements[CCCS]:
    def getIndexIVS(Vsource):
        for i in range(len(cirElements[IVS])):
            if cirElements[IVS][i].name == Vsource:
                return i
    if cccs.n1 != 'GND':
        M[numbered_nodes[cccs.n1]][num_nodes+getIndexIVS(cccs.vsource)]-=cccs.value
    if cccs.n2 != 'GND':
        M[numbered_nodes[cccs.n2]][num_nodes+getIndexIVS(cccs.vsource)]+=cccs.value

try:
    x = np.linalg.solve(M, b)
    cirCurrents = []
    cirNodes = []
    # Formatting Output Data
    #print("M:\n")
    #print(M)
    #print("b:\n")
    #print(b)
    #print("x:\n")
    #print(x)
    for n in nodes:
        cirNodes.append("Voltage at "+ str(n))
    for v in cirElements[IVS]:
        cirCurrents.append("Current in "+v.name)
    for v in cirElements[VCVS]:
        cirCurrents.append("Current in "+v.name)
    for v in cirElements[CCVS]:
        cirCurrents.append("Current in "+v.name)
    # Printing output in table format
    print(pd.DataFrame(x, cirNodes+cirCurrents, columns=['Amplitude']))

except np.linalg.LinAlgError:
    print("Matrix M is singular. Hence not solvable.\n")
    sys.exit()


