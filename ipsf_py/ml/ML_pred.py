import pandas as pd
import numpy as np
from numpy.random import seed
seed(1)
from sklearn.model_selection import train_test_split 
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from sklearn import svm
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error
from math import sqrt
from scipy.stats import pearsonr
from pickle import dump
from pickle import load

def GET_PI(list_ep):
   ndata = len(list_ep)
   sum1 = 0.0
   sum2 = 0.0
   for i in range(ndata):
      ei = list_ep[i][0]
      pi = list_ep[i][1]
      for j in range(i+1, ndata):
         ej = list_ep[j][0]
         pj = list_ep[j][1]
         wij = abs(ej - ei)
         #if abs(pj-pi) < 0.001:
         #   cij = 0.0
         #elif abs(ej-ei) < 0.001:
         #   cij = 0.0
         #elif (ej-ei)/(pj-pi) < 0.0:
         #   cij = -1.0
         #elif (ej-ei)/(pj-pi) > 0.0:
         #   cij = 1.0
         #else:
         #   print "CAUTION: situation not considered, ei=%g, pi=%g, ej=%g, pj=%g"%(ei,pi,ej,pj)
         #   raise SystemExit
         #sum1 += wij * cij
         sum2 += wij
         if (ej-ei)*(pj-pi) > 0.0:
            sum1 = sum1 + wij
         elif (ej-ei)*(pj-pi) < 0.0:
            sum1 = sum1 - wij
   pi = sum1 / sum2
   return pi


def ml_save(kres_path, expt_path, model_lis):
   df_res = pd.read_csv(kres_path, index_col = 'molecule')
   df_exp = pd.read_csv(expt_path, index_col = 'ID')

   for i in df_res.index:
      if i in df_exp.index:
         df_res.loc[i, 'exp'] = df_exp.loc[i, 'exp']

   df_combine = df_res.dropna()

   df_combine=df_combine[df_combine.exp < 0]

   # X = df_combine.drop(columns = ['exp']) # for normal python 
   X = df_combine.drop(columns = ['exp'])

   y = df_combine[['exp']]


   #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state = None)
   # Initiate training variable and target value variable
   X_train = X
   y_train = y


   # Initiate stat matrix
   rmse = []
   mae = []
   corr = []
   corr2 = []
   p_value = []
   pi = []

   act_model = {}

   for mol in model_lis:
      # Initiate ML model instance
      if mol == 1:
         reg = linear_model.LinearRegression()
         name = 'LR'
         act_model[mol] = name
      elif mol == 2:
         reg = linear_model.Lasso(alpha=0.1)
         name = 'LASSO'
         act_model[mol] = name
      elif mol == 3:
         reg = linear_model.BayesianRidge()
         name = 'Bayesian'
         act_model[mol] = name
      elif mol == 4:
         reg = svm.SVR()
         name = 'SVR'
         act_model[mol] = name
      elif mol == 5:
         reg = RandomForestRegressor(max_depth=5, random_state=4)
         name = 'RF'
         act_model[mol] = name
      elif mol == 6:
         reg = AdaBoostRegressor(random_state=4, n_estimators=100)
         name = 'Adaboost'
         act_model[mol] = name
      elif mol == 7:
         reg = GradientBoostingRegressor(random_state=4)
         name = 'GBDT'
         act_model[mol] = name
      elif mol == 8:
         reg = MLPRegressor(random_state=4, max_iter=500, hidden_layer_sizes = (64,64,32), activation = 'relu')
         name = 'MLP'
         act_model[mol] = name

      # Train model
      reg.fit(X_train,y_train)

      # Save ML model
      with open('%s.pkl'%name, 'wb') as f:
         dump(reg, f, protocol=5)

      # Calc model error
      test_pred = reg.predict(X_train)
      rmse.append(sqrt(mean_squared_error(y_train, test_pred)))
      mae.append(mean_absolute_error(y_train, test_pred))
      print('\n%s'%name)
      a, b = pearsonr(y_train['exp'].values, test_pred.ravel())
      corr.append(a)
      p_value.append(b)
      corr2.append(np.square(a))
      temp = list(zip(y_train['exp'].values, test_pred.ravel()))
      pi.append(GET_PI(temp))

   # Summarize trained model performance
   value = {'RMSE':rmse, 'MAE':mae, 'CORR':corr, 'CORR2':corr2, 'p-value':p_value, 'PI':pi}
   estimation=pd.DataFrame(data=value)
   file_name = 'evaluation_ML_{}.csv'.format('train_all')
   estimation.to_csv(file_name, index = False, header = True)
   with open('MODEL_LIST.csv', 'a') as f:
      for m in model_lis:
         f.write("%s:,%s\n"%(str(m),act_model[m]))


def ml_bs(kres_path, expt_path, model_lis, split_ratio, n_bs):
   df_res = pd.read_csv(kres_path, index_col = 'molecule')
   df_exp = pd.read_csv(expt_path, index_col = 'ID')

   for i in df_res.index:
      if i in df_exp.index:
         df_res.loc[i, 'exp'] = df_exp.loc[i, 'exp']

   df_combine = df_res.dropna()

   df_combine=df_combine[df_combine.exp < 0]

   X = df_combine.drop(columns = ['exp'])

   y = df_combine[['exp']]

   for n in range(1,n_bs+1):
      X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=split_ratio, random_state = None)

      # Initiate stat matrix
      rmse = []
      mae = []
      corr = []
      corr2 = []
      p_value = []
      pi = []

      act_model = {}

      for mol in model_lis:
         # Initiate ML model instance
         if mol == 1:
            reg = linear_model.LinearRegression()
            name = 'LR'
            act_model[mol] = name
         elif mol == 2:
            reg = linear_model.Lasso(alpha=0.1)
            name = 'LASSO'
            act_model[mol] = name
         elif mol == 3:
            reg = linear_model.BayesianRidge()
            name = 'Bayesian'
            act_model[mol] = name
         elif mol == 4:
            reg = svm.SVR()
            name = 'SVR'
            act_model[mol] = name
         elif mol == 5:
            reg = RandomForestRegressor(max_depth=5, random_state=4)
            name = 'RF'
            act_model[mol] = name
         elif mol == 6:
            reg = AdaBoostRegressor(random_state=4, n_estimators=100)
            name = 'Adaboost'
            act_model[mol] = name
         elif mol == 7:
            reg = GradientBoostingRegressor(random_state=4)
            name = 'GBDT'
            act_model[mol] = name
         elif mol == 8:
            reg = MLPRegressor(random_state=4, max_iter=500, hidden_layer_sizes = (64,64,32), activation = 'relu')
            name = 'MLP'
            act_model[mol] = name
         
         # Train model
         reg.fit(X_train,y_train)

         # Save ML model
         with open('%s.pkl'%name, 'wb') as f:
            dump(reg, f, protocol=5)

         # Calc model error
         test_pred = reg.predict(X_test)
         rmse.append(sqrt(mean_squared_error(y_test, test_pred)))
         mae.append(mean_absolute_error(y_test, test_pred))
         a, b = pearsonr(y_test['exp'].values, test_pred.ravel())
         corr.append(a)
         p_value.append(b)
         corr2.append(np.square(a))
         temp = list(zip(y_test['exp'].values, test_pred.ravel()))
         pi.append(GET_PI(temp))

      # Summarize trained model performance
      value = {'RMSE':rmse, 'MAE':mae, 'CORR':corr, 'CORR2':corr2, 'p-value':p_value, 'PI':pi}
      estimation=pd.DataFrame(data=value)
      file_name = 'evaluation_ML_{}.csv'.format(n)
      estimation.to_csv(file_name, index = False, header = True)

   with open('MODEL_LIST.csv', 'a') as f:
      for m in model_lis:
         f.write("%s:,%s"%(str(m),act_model[m]))


def ml_load(kres_path, model_lis, model_root_path, **kwargs):
   df_res = pd.read_csv(kres_path, index_col = 'molecule')
   expt_path = kwargs.get('expt_path', None)
   if expt_path is None:
      expt_flag = False
   else:
      expt_flag = True
      df_exp = pd.read_csv(expt_path, index_col = 'ID')
      for i in df_res.index:
         #root_i = "_".join(i.split("_")[:-1])
         if i in df_exp.index:
            df_res.loc[i, 'exp'] = df_exp.loc[i, 'exp']

   df_combine = df_res.dropna()

   if expt_flag:
      df_combine=df_combine[df_combine.exp < 0]
      X = df_combine.drop(columns = ['exp'])
      y = df_combine[['exp']]
   else:
      X = df_combine
      y = pd.DataFrame(index=df_combine.index)

   act_model = {}
   out = y

   for mol in model_lis:
      # Initiate ML model instance
      if mol == 1:
         with open('%s/LR.pkl'%model_root_path, 'rb') as f:
            reg = load(f)
         name = 'LR'
         act_model[mol] = name
      elif mol == 2:
         with open('%s/LASSO.pkl'%model_root_path, 'rb') as f:
            reg = load(f)
         name = 'LASSO'
         act_model[mol] = name
      elif mol == 3:
         with open('%s/Bayesian.pkl'%model_root_path, 'rb') as f:
            reg = load(f)
         name = 'Bayesian'
         act_model[mol] = name
      elif mol == 4:
         with open('%s/SVR.pkl'%model_root_path, 'rb') as f:
            reg = load(f)
         name = 'SVR'
         act_model[mol] = name
      elif mol == 5:
         with open('%s/RF.pkl'%model_root_path, 'rb') as f:
            reg = load(f)
         name = 'RF'
         act_model[mol] = name
      elif mol == 6:
         with open('%s/Adaboost.pkl'%model_root_path, 'rb') as f:
            reg = load(f)
         name = 'Adaboost'
         act_model[mol] = name
      elif mol == 7:
         with open('%s/GBDT.pkl'%model_root_path, 'rb') as f:
            reg = load(f)
         name = 'GBDT'
         act_model[mol] = name
      elif mol == 8:
         with open('%s/MLP.pkl'%model_root_path, 'rb') as f:
            reg = load(f)
         name = 'MLP'
         act_model[mol] = name

      test_pred = reg.predict(X)
      out[mol] = test_pred

   out.to_csv('PRED.csv', index = True, header= True)