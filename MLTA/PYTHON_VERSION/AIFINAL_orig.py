#Code by Filipe Ossege

#Imported packages
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # <--- This is important for 3d plotting
from scipy.interpolate import griddata
import json
# To install the scikit-learn (sklearn) package in the Spyder software, use the following commands
# first install the miniconda anaconda distribution in your machine
# Open the anaconda prompt from miniconda and create an environment with the command
# "conda create -n spyder-env -y"
# Activate the environment with
# "conda activate spyder-env"
# Then install the scikit-learn package
# "pip install scikit-learn"
# Then, open the Spyder software. Click in "Tools" => "Preferences" => "Python interpreter"
# Click on "use the following Python interpreter" and choose ""spyder-env/python.ex"
# Click on apply. Then, click in "restart kernel" in the options button on the right corner
# The package can now be used

# To install python in ubuntu, open the terminal and type
# "sudo apt install python3"
# To install the sklearn, matplotlib and scipy package, type
# "pip3 install sklearn matplotlib scipy"
# Then, to run the routine, type in terminal
# "python3 AIFINAL.py"

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor

# Dataset
A = np.loadtxt("database.dat") # Varibles in columns (phi, H_0 [A/m], f [Hz], a [m]) and  (Tc[°C], t[s])
X = A[:,0:4] #phi, H_0 [A/m], f [Hz], a [m]
Y = A[:,4:6] #Tc[°C], t[s]

# Division in test data
XTrain, XTest, YTrain, YTest = train_test_split(X, Y, test_size=0.2)
YTrain1 = YTrain[:,0]
YTrain2 = YTrain[:,1]
YTest1 = YTest[:,0]
YTest2 = YTest[:,1]

# Standardize
scaler = StandardScaler()
XTrain_scaled = scaler.fit_transform(XTrain)
XTest_scaled = scaler.transform(XTest)
scaler_Y1 = StandardScaler()
YTrain1_scaled = scaler_Y1.fit_transform(YTrain1.reshape(-1, 1)).ravel()
scaler_Y2 = StandardScaler()
YTrain2_scaled = scaler_Y2.fit_transform(YTrain2.reshape(-1, 1)).ravel()

# first model (temperature)
Mdl1 = MLPRegressor(hidden_layer_sizes=(161, 10, 300),
                    activation='logistic', #logistic = sigmoid
                    alpha=5.0754e-06,
                    max_iter=3000,
                    solver='lbfgs')
Mdl1.fit(XTrain_scaled, YTrain1_scaled)

# Second model (time)
Mdl2 = MLPRegressor(hidden_layer_sizes=(10, 293),
                    activation='logistic', #logistic = sigmoid
                    alpha=2.1568e-08,
                    max_iter=3000,
                    solver='lbfgs')
Mdl2.fit(XTrain_scaled, YTrain1_scaled)

testPredictions1 = scaler_Y1.inverse_transform(Mdl1.predict(XTest_scaled).reshape(-1, 1)).ravel()
testPredictions2 = scaler_Y2.inverse_transform(Mdl2.predict(XTest_scaled).reshape(-1, 1)).ravel()

plt.figure(figsize=(4, 4))
plt.plot(YTest1, testPredictions1, ".")
plt.plot(YTest1, YTest1)
plt.xlabel("True Temperature (°C)")
plt.ylabel("Predicted Temperature (°C)")
plt.savefig('Temperature.png', format='png', dpi=1200, bbox_inches='tight')

plt.figure(figsize=(4, 4))
plt.plot(YTest2, testPredictions2, ".")
plt.plot(YTest2, YTest2)
plt.xlabel("True Temperature (°C)")
plt.ylabel("Predicted Temperature (°C)")
plt.savefig('TIME.png', format='png', dpi=1200, bbox_inches='tight')

R1 = np.concatenate((XTest, YTest1.reshape(-1, 1), testPredictions1.reshape(-1, 1)), axis=1)
R2 = np.concatenate((XTest, YTest2.reshape(-1, 1), testPredictions2.reshape(-1, 1)), axis=1)

###############################################################################
'''Predict temperature and time for custom input'''
phi = {phi} # % 
H_0 = {H_0} #[A/m] 
f = {freq} # [Hz]
a = {rad} #[m]
# Input given by the user
user_input = np.array([[phi, H_0, f, a]])

# Maximum values of the model
X_max = np.max(XTrain, axis=0) 
# Mininum values of the model
X_min = np.min(XTrain, axis=0) 
# Check if the input is under the validity range of the model
variables = ["phi", "H_0", "f", "a"]
for i, var in enumerate(variables):
    if not (X_min[i] <= user_input[0][i] <= X_max[i]):
        print(f"Caution: {var} is out of the training range!")
        print(f" - Variable {var}: value {user_input[0][i]} out of range [{X_min[i]}, {X_max[i]}]")
        break
else:
    print("All values are within training range.")

    # Scaled to zero mean and unitary standart deviation
    user_input_scaled = scaler.transform(user_input)
    # Predicts the value
    predicted_temp_scaled = Mdl1.predict(user_input_scaled)
    predicted_time_scaled = Mdl2.predict(user_input_scaled)
    # Transforms back to the original scale
    predicted_temp = scaler_Y1.inverse_transform(predicted_temp_scaled.reshape(-1, 1)).ravel()
    predicted_time = scaler_Y2.inverse_transform(predicted_time_scaled.reshape(-1, 1)).ravel()
    print(f"Predicted Temperature: {predicted_temp[0]:.2f} °C")
    print(f"Predicted Time: {predicted_time[0]:.2f} s")
    
resultado = {
    "predicted_temp": float(predicted_temp[0]),
    "predicted_time": float(predicted_time[0])
}

with open("resultado.json", "w") as f:
    json.dump(resultado, f)
    
with open("resultado.json", "r") as f:
    data = json.load(f)

predicted_temp = data["predicted_temp"]
predicted_time = data["predicted_time"]