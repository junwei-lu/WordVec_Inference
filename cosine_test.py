import h5py
import torch

diag = torch.diagonal

def cosine_test_diag(i, Pi, diagP, diagPDP, PDPi, P1, p_w, myAi, n ,T1):
    """
Cosine_test_diag = function(testid = 1,P,P1,n,T1,p_w){
  # testid denote the row of P which will be tested by estimating it's covariance
  # P is defined in the paper around Formula (2.7) as U%*%t(U)
  # P1 is the rowsum of P
  # n is the number of patients
  # T1 is defined in the paper around  Formula:(3.9) as T1=T*q-q*(q+1)/2
  # p_w is defined in the paper in Algorithm 1
}
    """
    P1i = P1[i]
    p_wi = p_w[i]
    d = p_w.shape[0]
    a = Pi ** 2 / p_w
    f = Pi ** 2 / p_w**2
    b = Pi * (P1i - Pi) / p_w
    # Estimating Covariance of the testid's row
    # step1,2,3,4 denote the (1) (2) (3) (4) part of covariance in the paper (3.11)&(3.12)

    # step1
    EstCov = P1**2 * (2*p_wi-1)
    tmp = Pi*P1
    EstCov = EstCov - tmp*2 + diagPDP * (1-p_wi) + (Pi**2)/p_wi
    EstCov = EstCov / p_wi
    # step 2
    repdiagPDP = (-diagPDP[i]).clone().detach().repeat(d).to("cuda:0")
    step2 = repdiagPDP - a*2 + 1/p_w*diagPDP[i] + f 
    EstCov = EstCov + step2
    # step 3
    step3 = (2*(P1i**2)).clone().detach().repeat(d).to("cuda:0")
    step3 = step3 - b*2 - (1/p_w*P1i**2) ## changed here
    EstCov = EstCov + step3
    # step 4
    step4 = P1*4*P1i
    step4 = step4 - PDPi*2
    step4[i] = step4[i] - P1i**2/p_wi*2
    tmp = diagP/p_w
    step4 = step4 - tmp*2*P1i
    step4[i] = step4[i] + diagPDP[i]/p_wi*2
    EstCov = EstCov + step4
    # merge
    EstCov = EstCov/(2*n*T1)
    # PMI to Cos
    EstCov = EstCov*(myAi**2)
    return EstCov

print("Loading P")
fn = "P.hdf5"
with h5py.File(fn, "r") as f:
    P = torch.tensor(f["pmi"][...])
print("Loading PDP")
fn = "PDP.hdf5"
with h5py.File(fn, "r") as f:
    PDP = torch.tensor(f["pmi"][...])
print("Loading myA")
fn = "myA.hdf5"
with h5py.File(fn, "r") as f:
    myA = torch.tensor(f["pmi"][...])
print("Loading p_w")
fn = "p_w.hdf5"
with h5py.File(fn, "r") as f:
    p_w = torch.tensor(f["pmi"][...]).t()[:,0]

d = p_w.shape[0]
device = torch.device("cuda:0")
n = 12800000 # input the n computed in the pre.R
T1 = 3880903. # input the T1 computed in the pre.R
diagP = diag(P).to(device)
diagPDP = diag(PDP).to(device)
P1 = torch.sum(P,1).to(device)
p_w = p_w.to(device)
cosineTest = torch.zeros(d,d)
print("Calling cosine_test_diag")

for i in range(d):
   Pi = P[:,i].to(device)
   PDPi = PDP[i,].to(device)
   myAi = myA[i,].to(device)
   x = cosine_test_diag(i, Pi, diagP, diagPDP, PDPi, P1, p_w, myAi, n ,T1) 
   cosineTest[i,] = x
print("Done")
fn = "cosineTest.hdf5"
with h5py.File(fn, "w") as f:
    f.create_dataset("cosineTest",data=cosineTest)