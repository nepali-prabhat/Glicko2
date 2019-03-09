import math
r = 1500
RD = 200
sigma = 0.06
tau = 0.5
opponents={"r":[1400,1550,1700], "RD":[30,100,300]}
s = [1,0,0]

def glickon_meu(r):
    return (r-1500)/173.7178

def glickon_phi(rd):
    return rd/173.7178

def g(phi):
    value = 1+((3*(phi**2))/math.pi**2)
    return math.pow(value,-1/2)
def E(meu,meu_j,phi_j):
    value = 1+(math.exp(-g(phi_j)*(meu-meu_j)))
    return math.pow(value,-1)

meu = glickon_meu(r)
print('meu=', meu)
phi = glickon_phi(RD)
print('phi=', phi)

meu_opponents = [glickon_meu(x) for x in opponents['r']]

phi_opponents = [glickon_phi(x) for x in opponents['RD']]

def sum_value_v(meu,meu_j,phi_j):
    e = E(meu,meu_j,phi_j) 
    return g(phi_j)**2*e*(1-e)

def v(meu,meu_opponents,phi_opponents):
    value = 0;
    for j in range(len(meu_opponents)):
        value += sum_value_v(meu, meu_opponents[j], phi_opponents[j])
    return math.pow(value,-1)

value_of_v = v(meu,meu_opponents,phi_opponents)
print("v=",value_of_v)

def sum_value_delta(s_j,meu,meu_j,phi_j):
    return g(phi_j)*(s_j - E(meu,meu_j,phi_j))

def delta(v,meu,meu_opponents,phi_opponents,s):
    value = 0;
    for j in range(len(meu_opponents)):
        value += sum_value_delta(s[j],meu,meu_opponents[j],phi_opponents[j])
    value=v*value
    return value

value_of_delta = delta(value_of_v,meu,meu_opponents,phi_opponents,s)

print("delta=",value_of_delta)

a = math.log(sigma**2)

f = lambda x:( ( math.exp(x)*(value_of_delta**2 - phi**2 - value_of_v - math.exp(x)) ) / (2* (phi**2 + value_of_v + math.exp(x))**2 ) ) - ((x-a)/tau**2)

e = 0.000001

def new_sigma(a, value_of_delta,value_of_v,phi,tau):
    A = a
    if value_of_delta**2 > (phi**2 + value_of_v):
        B = math.log(value_of_delta**2 - phi**2- value_of_v)
    else:
        k = 1;
        while f(a-k*math.sqrt(tau**2))<0:
               k+=1;
        B =  a-k*math.sqrt(tau**2)
    f_a = f(A)
    f_b = f(B)
    while math.fabs(B-A) > e:
        C = A + (A-B)* f_a/(f_b-f_a)
        f_c = f(C)
        if f_c*f_b < 0:
            A = B
            f_a = f_b
        else:
            f_a = f_a/2
        B = C 
        f_b = f_c
    sigma_prime = math.exp(A/2)
    return sigma_prime
sigma_prime = new_sigma(a, value_of_delta,value_of_v,phi,tau)
print("sigma prime=",sigma_prime)
def phi_star():
    return math.pow(phi**2+sigma_prime**2,1/2)
print("phi star=",phi_star())

def new_phi(phi, new_sigma, value_of_v):
    val =  math.pow(phi_star(),-1) + math.pow(value_of_v,-1)
    return math.pow(val, -1/2)

phi_prime = new_phi(phi, new_sigma, value_of_v)
print("phi prime=",phi_prime)
def sum_value_new_meu(phi_j,meu,meu_j,s_j):

    return g(phi_j)*(s_j-E(meu,meu_j,phi_j))

def new_meu(meu, new_phi, phi_opponents, meu_opponents, s):
    value = 0
    for j in range(len(phi_opponents)):
        value += sum_value_new_meu(phi_opponents[j],meu,meu_opponents[j],s[j])
    value = phi_prime**2 * value
    value = value + meu
    return value
meu_prime = new_meu(meu, new_phi, phi_opponents, meu_opponents, s)
print("meu prime=",meu_prime)

def r_from_glicko_scale(meu_prime):
    return meu_prime*173.7178 + 1500

def RD_from_glicko_scale(phi_prime):
    return 173.7178*phi_prime

def sigma_from_glicko_scale(sigma_prime):
    return sigma_prime

print("\n\nnew values: ")
print("r=",r_from_glicko_scale(meu_prime))
print("RD=",RD_from_glicko_scale(phi_prime))
print("sigma=",sigma_from_glicko_scale(sigma_prime))