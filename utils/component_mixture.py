import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

# Константы:
R0 = 8.31


def initialize_mixture_component_dispersed(mixture_reference, name, T_base, mu, rho, dH0, eps_dk, sigma, T_grid, Cp_grid, H_grid):
    if name in mixture_reference.keys():
        return Component_Dispersed(name, T_base, mu, rho, dH0, eps_dk, sigma, T_grid, Cp_grid, H_grid)
    else:
        print('mixture dispersed component with name ', name, ' not defined')
        return None


def initialize_mixture_component_gas(mixture_reference, name, T_base, mu, rho, dH0, eps_dk, sigma, T_grid, Cp_grid, H_grid, viscosity_grid, heat_conductivity_grid, diffusivity_grid):
    if name in mixture_reference.keys():
        return Component_Gas(name, T_base, mu, rho, dH0, eps_dk, sigma, T_grid, Cp_grid, H_grid, viscosity_grid, heat_conductivity_grid, diffusivity_grid)
    else:
        print('mixture gas component with name ', name, ' not defined')
        return None


class Component:
    def __init__(self, name, T_base, mu, rho, dH0, eps_dk, sigma, T_grid, Cp_grid, H_grid):
        self.name = name
        self.mu = mu
        self.rho = rho
        self.dH0 = dH0
        self.eps_dk = eps_dk
        self.sigma = sigma
        self.T_base = T_base
        self.T_grid = T_grid
        self.Cp_grid = Cp_grid
        self.H_grid = H_grid

    def interpolated_value(self, x, x0, x1, y0, y1):
        return y0 + (y1 - y0) / (x1 - x0) * (x - x0)

    def get_value_from_grid(self, T, grid):
        if T < self.T_base:
            # in case zero included in T_grid (I dont remember)
            it = np.searchsorted(self.T_grid, T, side='left')
            return grid[it - 1]
        if T > max(self.T_grid):
            return grid[len(grid) - 1]
        else:
            it = np.searchsorted(self.T_grid, T, side='left')
            T0 = self.T_grid[it - 1]
            T1 = self.T_grid[it]
            y_0 = grid[it - 1]
            y_1 = grid[it]
            return self.interpolated_value(T, T0, T1, y_0, y_1)

    def get_Cp(self, T):
        return self.get_value_from_grid(T, self.Cp_grid)

    # def get_Cp(self, T):
    #     if T < min(self.T_grid):
    #         return self.Cp_grid[0]
    #     if T > max(self.T_grid):
    #         return self.Cp_grid[len(self.Cp_grid)-1]
    #     else:
    #         it = np.searchsorted(self.T_grid, T, side='left')
    #         T0 = self.T_grid[it-1]
    #         T1 = self.T_grid[it]
    #         Cp_0 = self.Cp_grid[it-1]
    #         Cp_1 = self.Cp_grid[it]
    #         return self.interpolated_value(T, T0, T1, Cp_0, Cp_1)

    def get_H(self, T):
        # print('inside get_H for mixture component')
        if T < self.T_grid[0]:
            # print('inside less loop')
            it = np.searchsorted(self.T_grid, T, side='left')
            return self.Cp_grid[it - 1]*T
            # return self.Cp_grid[0]*T
        max_index = len(self.T_grid) - 1
        if T > self.T_grid[max_index]:
            # print('inside more loop')
            return self.H_grid[max_index] + self.Cp_grid[max_index] * (T - self.T_grid[max_index])
        else:
            it = np.searchsorted(self.T_grid, T, side='left')
            T0 = self.T_grid[it-1]
            T1 = self.T_grid[it]
            H_0 = self.H_grid[it-1]
            H_1 = self.H_grid[it]
            # print('inside general loop', it, T0, T1, H_0, H_1)
            return self.interpolated_value(T, T0, T1, H_0, H_1)

    def get_mu(self):
        return self.mu

    def get_dH0(self):
        return self.dH0


class Component_Dispersed(Component):
    def __init__(self, name, T_base, mu, rho, dH0, eps_dk, sigma, T_grid, Cp_grid, H_grid):
        super().__init__(name, T_base, mu, rho, dH0, eps_dk, sigma, T_grid, Cp_grid, H_grid)
        self.columns_mix = ['T [K]', 'Cp [Дж/кг-К]', 'H [Дж/кг]']
        self.data = pd.DataFrame(data=np.array([self.T_grid, self.Cp_grid, self.H_grid]).transpose(), index=self.T_grid, columns=self.columns_mix)
        # self.data.drop(index=0, inplace=True)
        print('mixture dispersed component ', name, ' is initialized')

    def get_rho(self):
        return self.rho


class Component_Gas(Component):
    def __init__(self, name, T_base, mu, rho, dH0, eps_dk, sigma, T_grid, Cp_grid, H_grid, viscosity_grid, heat_conductivity_grid, diffusivity_grid):
        super().__init__(name, T_base, mu, rho, dH0, eps_dk, sigma, T_grid, Cp_grid, H_grid)
        self.viscosity_grid = viscosity_grid
        self.heat_conductivity_grid = heat_conductivity_grid
        self.diffusivity_grid = diffusivity_grid
        self.columns_mix = ['T [K]', 'Cp [Дж/кг-К]', 'H [Дж/кг]', 'Mu_visc [kg/m-s]', 'Lambda [W/m-K]', 'D [m^2/s]']
        self.data = pd.DataFrame(data=np.array([self.T_grid, self.Cp_grid, self.H_grid, self.viscosity_grid, self.heat_conductivity_grid, self.diffusivity_grid]).transpose(), index=self.T_grid, columns=self.columns_mix)
        # self.data.drop(index=0, inplace=True)
        print('mixture gas component ', name, ' is initialized')


    def get_viscosity(self, T):
        return self.get_value_from_grid(T, self.viscosity_grid)

    def get_heat_conductivity(self, T):
        return self.get_value_from_grid(T, self.heat_conductivity_grid)

    def get_diffusivity(self, T):
        return self.get_value_from_grid(T, self.diffusivity_grid)
