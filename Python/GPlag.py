import math
import torch
import gpytorch
from torch.nn.functional import softplus
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from gpytorch.constraints import Interval
from torch.distributions.multivariate_normal import MultivariateNormal
# set to training mode and train

# Define the custom kernel
class MTSRBF(gpytorch.kernels.Kernel):

    def __init__(self, sigma=1, alpha = 1, length = 1, timelag_init_value = 0, sigma_constraint = None, alpha_constraint = None, length_constraint = None, timelag_lower_bound = None, timelag_upper_bound = None):
        super(MTSRBF, self).__init__()
        self.raw_alpha = torch.nn.Parameter(torch.tensor(alpha, dtype=torch.float32).log())
        self.register_parameter(name="raw_alpha", parameter=self.raw_alpha)
        if alpha_constraint is None:
            alpha_constraint = gpytorch.constraints.Positive()
        # register the constraint
        self.register_constraint("raw_alpha", alpha_constraint)

        self.raw_sigma = torch.nn.Parameter(torch.tensor(sigma, dtype=torch.float32).log())
        self.register_parameter(name="raw_sigma", parameter=self.raw_sigma)

        if sigma_constraint is None:
            sigma_constraint = gpytorch.constraints.Positive()
        # register the constraint
        self.register_constraint("raw_sigma", sigma_constraint)

        
        self.raw_length = torch.nn.Parameter(torch.tensor(length, dtype=torch.float32).log())
        self.register_parameter(name="raw_length", parameter=self.raw_length)

        # set the parameter constraint to be positive, when nothing is specified
        if length_constraint is None:
            length_constraint = gpytorch.constraints.Positive()
        
        # register the constraint
        self.register_constraint("raw_length", length_constraint)

        # Register the timelag parameter without a transformation
        self.register_parameter(name="raw_timelag", parameter=torch.nn.Parameter(torch.tensor(timelag_init_value, dtype=torch.float32)))

        # Create a custom constraint for the timelag parameter
        if timelag_lower_bound is not None and timelag_upper_bound is not None:
            self.register_constraint("raw_timelag", Interval(lower_bound=torch.tensor(timelag_lower_bound, dtype=torch.float32), upper_bound=torch.tensor(timelag_upper_bound, dtype=torch.float32)))

    @property
    def alpha(self):
        return self.raw_alpha_constraint.transform(self.raw_alpha)

    @alpha.setter
    def alpha(self, value):
        return self._set_alpha(value)

    def _set_alpha(self, value):
        if not torch.is_tensor(value):
            value = torch.as_tensor(value).to(self.raw_alpha)
        # when setting the paramater, transform the actual value to a raw one by applying the inverse transform
        self.initialize(raw_alpha=self.raw_alpha_constraint.inverse_transform(value))

    @property
    def sigma(self):
        return self.raw_sigma_constraint.transform(self.raw_sigma)

    @sigma.setter
    def sigma(self, value):
        return self._set_sigma(value)

    def _set_sigma(self, value):
        if not torch.is_tensor(value):
            value = torch.as_tensor(value).to(self.raw_sigma)
        # when setting the paramater, transform the actual value to a raw one by applying the inverse transform
        self.initialize(raw_sigma=self.raw_sigma_constraint.inverse_transform(value))

    @property
    def length(self):
        # when accessing the parameter, apply the constraint transform
        return self.raw_length_constraint.transform(self.raw_length)

    @length.setter
    def length(self, value):
        return self._set_length(value)

    def _set_length(self, value):
        if not torch.is_tensor(value):
            value = torch.as_tensor(value).to(self.raw_length)
        # when setting the paramater, transform the actual value to a raw one by applying the inverse transform
        self.initialize(raw_length=self.raw_length_constraint.inverse_transform(value))

    def forward(self, x1, x2, n1 = None, n2 = None, delta = None, **params):
        ngroup_x2 = x2.shape[0]-x1.shape[0]
        if ngroup_x2 == 0:
            x1_ = x1 + torch.concat([torch.zeros(n1,1),self.raw_timelag * torch.ones(n2,1)]) 
            x2_ = x2 + torch.concat([torch.zeros(n1,1),self.raw_timelag * torch.ones(n2,1)])
            if delta is None:
                blk = np.ones((n1, n2), int)
                delta = torch.tensor(np.block([
                [np.zeros((n1, n1)), blk],
                [blk, np.zeros((n2, n2))]
            ]), dtype=torch.float32)
        else:
            ntest = x2.shape[0] // 2
            x1_ = x1 + torch.concat([torch.zeros(ntest,1),self.raw_timelag * torch.ones(ntest,1)]) 
            x2_ = x2 + torch.concat([torch.zeros(n1,1),self.raw_timelag * torch.ones(n2,1),torch.zeros(ntest,1),self.raw_timelag * torch.ones(ntest,1)])
            if delta is None:
                blk_1 = np.ones((ntest, n2), int)
                blk_2 = np.ones((ntest, ntest), int)
                blk_3 = np.ones((ntest, n1), int)
                delta = torch.tensor(np.block([
                [np.zeros((ntest, n1)), blk_1, np.zeros((ntest, ntest)), blk_2],
                [blk_3, np.zeros((ntest, n2)), blk_2, np.zeros((ntest, ntest))]
            ]), dtype=torch.float32)

        Dnom = self.covar_dist(x1_, x2_, square_dist=True)
        A = (self.alpha)*(delta)+1
        D = Dnom.div(-A)
        C = torch.exp(D * self.length)
        jitter = torch.tensor(1e-4)
        cov_mat = self.sigma * C.div(torch.sqrt(A)) + torch.eye(n1+n2) * jitter
        return cov_mat

class MTSE(gpytorch.kernels.Kernel):

    def __init__(self, sigma=1, alpha = 1, length = 1, timelag_init_value = 0, sigma_constraint = None, alpha_constraint = None, length_constraint = None, timelag_lower_bound = None, timelag_upper_bound = None):
        super(MTSE, self).__init__()
        self.raw_alpha = torch.nn.Parameter(torch.tensor(alpha, dtype=torch.float32).log())
        self.register_parameter(name="raw_alpha", parameter=self.raw_alpha)
        if alpha_constraint is None:
            alpha_constraint = gpytorch.constraints.Positive()
        # register the constraint
        self.register_constraint("raw_alpha", alpha_constraint)

        self.raw_sigma = torch.nn.Parameter(torch.tensor(sigma, dtype=torch.float32).log())
        self.register_parameter(name="raw_sigma", parameter=self.raw_sigma)

        if sigma_constraint is None:
            sigma_constraint = gpytorch.constraints.Positive()
        # register the constraint
        self.register_constraint("raw_sigma", sigma_constraint)

        
        self.raw_length = torch.nn.Parameter(torch.tensor(length, dtype=torch.float32).log())
        self.register_parameter(name="raw_length", parameter=self.raw_length)

        # set the parameter constraint to be positive, when nothing is specified
        if length_constraint is None:
            length_constraint = gpytorch.constraints.Positive()
        
        # register the constraint
        self.register_constraint("raw_length", length_constraint)

        # Register the timelag parameter without a transformation
        self.register_parameter(name="raw_timelag", parameter=torch.nn.Parameter(torch.tensor(timelag_init_value, dtype=torch.float32)))

        # Create a custom constraint for the timelag parameter
        if timelag_lower_bound is not None and timelag_upper_bound is not None:
            self.register_constraint("raw_timelag", Interval(lower_bound=torch.tensor(timelag_lower_bound, dtype=torch.float32), upper_bound=torch.tensor(timelag_upper_bound, dtype=torch.float32)))

    @property
    def alpha(self):
        return self.raw_alpha_constraint.transform(self.raw_alpha)

    @alpha.setter
    def alpha(self, value):
        return self._set_alpha(value)

    def _set_alpha(self, value):
        if not torch.is_tensor(value):
            value = torch.as_tensor(value).to(self.raw_alpha)
        # when setting the paramater, transform the actual value to a raw one by applying the inverse transform
        self.initialize(raw_alpha=self.raw_alpha_constraint.inverse_transform(value))

    @property
    def sigma(self):
        return self.raw_sigma_constraint.transform(self.raw_sigma)

    @sigma.setter
    def sigma(self, value):
        return self._set_sigma(value)

    def _set_sigma(self, value):
        if not torch.is_tensor(value):
            value = torch.as_tensor(value).to(self.raw_sigma)
        # when setting the paramater, transform the actual value to a raw one by applying the inverse transform
        self.initialize(raw_sigma=self.raw_sigma_constraint.inverse_transform(value))

    @property
    def length(self):
        # when accessing the parameter, apply the constraint transform
        return self.raw_length_constraint.transform(self.raw_length)

    @length.setter
    def length(self, value):
        return self._set_length(value)

    def _set_length(self, value):
        if not torch.is_tensor(value):
            value = torch.as_tensor(value).to(self.raw_length)
        # when setting the paramater, transform the actual value to a raw one by applying the inverse transform
        self.initialize(raw_length=self.raw_length_constraint.inverse_transform(value))

    def __call__(self, x1, x2=None, n1=None, n2=None, delta=None, diag=False, last_dim_is_batch=False, **params):
        # Call the parent class's __call__ method
        super().__call__(x1, x2, diag=diag, last_dim_is_batch=last_dim_is_batch, **params)

        # Pass n1, n2, and delta to the base kernel's forward method
        return self.forward(x1, x2, n1, n2, delta, **params)
        
    def forward(self, x1, x2, n1 = None, n2 = None, delta = None, **params):
        ngroup_x2 = x2.shape[0]-x1.shape[0]
        if ngroup_x2 == 0:
            x1_ = x1 + torch.concat([torch.zeros(n1,1),self.raw_timelag * torch.ones(n2,1)]) 
            x2_ = x2 + torch.concat([torch.zeros(n1,1),self.raw_timelag * torch.ones(n2,1)])
            if delta is None:
                blk = np.ones((n1, n2), int)
                delta = torch.tensor(np.block([
                [np.zeros((n1, n1)), blk],
                [blk, np.zeros((n2, n2))]
            ]), dtype=torch.float32)
        else:
            ntest = x2.shape[0] // 2
            x1_ = x1 + torch.concat([torch.zeros(ntest,1),self.raw_timelag * torch.ones(ntest,1)]) 
            x2_ = x2 + torch.concat([torch.zeros(n1,1),self.raw_timelag * torch.ones(n2,1),torch.zeros(ntest,1),self.raw_timelag * torch.ones(ntest,1)])
            if delta is None:
                blk_1 = np.ones((ntest, n2), int)
                blk_2 = np.ones((ntest, ntest), int)
                blk_3 = np.ones((ntest, n1), int)
                delta = torch.tensor(np.block([
                [np.zeros((ntest, n1)), blk_1, np.zeros((ntest, ntest)), blk_2],
                [blk_3, np.zeros((ntest, n2)), blk_2, np.zeros((ntest, ntest))]
            ]), dtype=torch.float32)

        Dnom = self.covar_dist(x1_, x2_, square_dist=False)
        A = (self.alpha)*(delta)+1
        C = torch.exp(-Dnom * self.length)
        jitter = torch.tensor(1e-4)
        cov_mat = self.sigma * C.div(A) + torch.eye(n1+n2) * jitter
        return cov_mat

class CustomScaleKernel(gpytorch.kernels.ScaleKernel):
    def __call__(self, x1, x2=None, n1=None, n2=None, delta=None, diag=False, last_dim_is_batch=False, **params):
        # Call the parent class's __call__ method
        super().__call__(x1, x2, diag=diag, last_dim_is_batch=last_dim_is_batch, **params)

        # Pass n1, n2, and delta to the base kernel's forward method
        return self.forward(x1, x2, n1, n2, delta, **params)

    def forward(self, x1, x2, n1, n2, delta, **params):
        # Call the base kernel's forward method with n1, n2, and delta
        return self.base_kernel.forward(x1, x2, n1, n2, delta, **params)

class GPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, kernel):
        super(GPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = CustomScaleKernel(kernel)

    def forward(self, x, n1, n2):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x, x, n1, n2)  # pass n1 and n2 here
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)

class MTSGP(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, kernel):
        super(MTSGP, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = kernel

    def forward(self, x, n1, n2):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x, x, n1, n2)  # pass n1 and n2 here
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)

# Wrap training, prediction and plotting from the ExactGP-Tutorial into a function,
# so that we do not have to repeat the code later on
def train_adam(model, likelihood, train_x, train_y, n1, n2, lr = 0.2, training_iter=200):
    # Use the adam optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)  # Includes GaussianLikelihood parameters
    # "Loss" for GPs - the marginal log likelihood
    mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)
    for i in range(training_iter):
        # Zero gradients from previous iteration
        optimizer.zero_grad()
        # Output from model
        output = model(train_x, n1 = n1, n2 = n2)
        # Calc loss and backprop gradients
        loss = -mll(output, train_y)
        loss.backward()
        optimizer.step()
        # print('Iter %d/%d - Loss: %.3f' % (i + 1, training_iter, loss.item()))
        # print(f"----- Iteration {i} -----")
        # for name, param in model.named_parameters():
        #     print(f"Parameter {name}: {param.item()}")