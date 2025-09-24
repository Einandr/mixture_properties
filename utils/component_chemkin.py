import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mendeleev import element
import re
import math

# Константы:
R0 = 8.314       # [Дж / моль-К] универсальная газовая постоянная
cal = 4184       # [кДж] термическая калория, используемая в Ansys Chemkin

# T_base = 100     # [K] температура для вычисления Ср с целью сведения в единое значение полной энтальпии при расчете через коэффициенты и через интеграл

# def initialize_mixture_component_dispersed(mixture_reference, name, T_base, mu, rho, dH0, T_grid, Cp_grid, H_grid):
#     if name in mixture_reference.keys():
#         return Component_Dispersed(name, T_base, mu, rho, dH0, T_grid, Cp_grid, H_grid)
#     else:
#         print('mixture dispersed component with name ', name, ' not defined')
#         return None
#
#
# def initialize_mixture_component_gas(mixture_reference, name, T_base, mu, rho, dH0, T_grid, Cp_grid, H_grid, viscosity_grid, heat_conductivity_grid, diffusivity_grid):
#     if name in mixture_reference.keys():
#         return Component_Gas(name, T_base, mu, rho, dH0, T_grid, Cp_grid, H_grid, viscosity_grid, heat_conductivity_grid, diffusivity_grid)
#     else:
#         print('mixture gas component with name ', name, ' not defined')
#         return None


class Component:
    def __init__(self, name, date, formula, phase, T_low, T_mid, T_high, atomic, a1_low, a2_low, a3_low, a4_low, a5_low, a6_low, a7_low, a1_high, a2_high, a3_high, a4_high, a5_high, a6_high, a7_high, dT, T_base, T0):
        self.name = name
        self.date = date
        self.formula = formula
        self.phase = phase
        self.T_low = T_low
        self.T_mid = T_mid
        self.T_high = T_high
        self.atomic = atomic
        self.T_range = np.array([T_low, T_mid, T_high])
        self.a = np.array([[a1_low, a2_low, a3_low, a4_low, a5_low, a6_low, a7_low], [a1_high, a2_high, a3_high, a4_high, a5_high, a6_high, a7_high]])
        self.dT = dT                                                        # интервал с которым будем создавать сетку значений по температуре
        self.T_base = T_base                                                # [K] температура для вычисления крайнего нижнего значения Ср с целью сведения в единое значение полной энтальпии при расчете через коэффициенты и через интеграл
        self.T0 = T0                                                        # должна быть 298.15 K - температура для вычисления энтальпии образования

        self.sp_list = (re.sub(r'(\d+[.,]?)', r'\1 ', self.formula))        # выделяем границу молей предыдущего компонента и имени следующего
        self.sp_list = (re.sub(r'([a-zA-Z]+)', r'\1 ', self.sp_list))       # выделяем границу компонента и его количества молей
        self.sp_list = re.split(r'\s', self.sp_list)                        # разбиваем строку по пробелам
        self.sp_list = [x for x in self.sp_list if x]

        self.sp_components = [x for ind, x in enumerate(self.sp_list) if ((ind % 2) != 1)]
        self.sp_moles = [x for ind, x in enumerate(self.sp_list) if ((ind % 2) != 0)]
        self.mu = 0                                                         # [г/моль]
        for ind, x in enumerate(self.sp_components):
            self.el = element(x)
            self.moles = float(self.sp_moles[ind])
            self.mu += self.el.atomic_weight * self.moles
        self.R = R0 / self.mu * 1000                                        # [Дж/(кг-К)]

        self.T_grid = np.arange(self.T_low, self.T_high + self.dT, self.dT)
        self.update_T_grid()

        self.vector_Cp = np.vectorize(self.Cp, excluded=['units'])
        self.Cp_grid = self.vector_Cp(T=self.T_grid, units='J_kg_K')




        # self.rho = rho
        # self.dH0 = dH0
        # self.T_base = T_base
        # self.T_grid = T_grid
        # self.Cp_grid = Cp_grid
        # self.H_grid = H_grid

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


    # YA COPIED NOW
    def _Cp_kcal_mole(self, T, a1, a2, a3, a4, a5):
        return a1 + a2*T + a3*pow(T, 2) + a4*pow(T, 3) + a5*pow(T, 4)

    def _Cp_J_kg_K(self, T, a1, a2, a3, a4, a5):
        # return self._Cp_kcal_mole(T, a1, a2, a3, a4, a5) * cal * 1000 / self.mu
        return self._Cp_kcal_mole(T, a1, a2, a3, a4, a5) * self.R

    def _H_kcal_mole(self, T, a1, a2, a3, a4, a5, a6):
        return a1*T + a2/2*pow(T, 2) + a3/3*pow(T, 3) + a4/4*pow(T, 4) + a5/5*pow(T, 5) + a6

    def _H_J(self, T, a1, a2, a3, a4, a5, a6):
        return self._H_kcal_mole(T, a1, a2, a3, a4, a5, a6) * self.R

    # Для вычисления H0_298 Нужно использовать только нижний температурный диапазон - потому как привязка а6 коэффициента идет к температурному диапазону
    # и он предназначен для вычисления просто энтальпии (не образования) в выбранном температурном диапазоне. Для нижнего температурного диапазона а6 будет равен
    # искусственной величине H0_0, то есть энтальпии образования при 0К, для верхнего диапазона для него вообще еще больше теряется физический смысл -
    # это как будто бы энтальпия образования при 0К, если бы не было нижнего температурного диапазона.
    # О таком вычислении H0_298 прямо говорится в Chemkin Input Manual пункт 2.3

    def H0_298_kcal_mole(self):
        return self._H_kcal_mole(self.T0, self.a[0, 0], self.a[0, 1], self.a[0, 2], self.a[0, 3], self.a[0, 4], self.a[0, 5])

    def H0_298_J(self):
        return self.H0_298_kcal_mole() * self.R

    def Cp(self, T, units):         # Эта функция используется для заполнения массива по запрашиваемой температуре
        # TODO добавить вычисление для Т_base 100К чтобы свести к одинаковому значению полные энтальпии при вычислении через коэффициенты и через интеграл Ср
        i = 0
        n_max = len(self.T_range) - 1
        n = 0
        while n < n_max:
            if T > self.T_range[n]:
                n += 1
            else:
                break
            i = n - 1
        if T >= self.T_range[n_max]:
            i = n_max - 1
            T = self.T_range[n_max]
        elif T <= self.T_range[0]:
            T = self.T_range[0]
        if units == 'kcal_mole':
            return self._Cp_kcal_mole(T, self.a[i, 0], self.a[i, 1], self.a[i, 2], self.a[i, 3], self.a[i, 4])
        elif units == 'J_kg_K':
            return self._Cp_J_kg_K(T, self.a[i, 0], self.a[i, 1], self.a[i, 2], self.a[i, 3], self.a[i, 4])
        else:
            raise ValueError('неверные единицы измерения Ср')

    def update_T_grid(self):
        # TODO добавить Т_base 100К
        self.T_grid = np.append(self.T_grid, self.T_mid)
        self.T_grid = np.unique(self.T_grid)
        return None










    # def H_From_Terra(self, T):
    #     i = 0
    #     n_max = len(self.T_range) - 1
    #     n = 0
    #     while n < n_max:
    #         if T > self.T_range[n]:
    #             n += 1
    #         else:
    #             break
    #         i = n - 1
    #     if T >= self.T_range[n_max]:
    #         i = n_max - 1
    #         T = self.T_range[n_max]
    #     elif T <= self.T_range[0]:
    #         T = self.T_range[0]
    #     # print(i)
    #     x = 0.0001*T
    #     return self.H(x, self.f_ranges[i, 1], self.f_ranges[i, 2], self.f_ranges[i, 3], self.f_ranges[i, 4], self.f_ranges[i, 5], self.f_ranges[i, 6])





    # def get_Cp(self, T):
    #     return self.get_value_from_grid(T, self.Cp_grid)
    #
    #
    # def get_H(self, T):
    #     # print('inside get_H for mixture component')
    #     if T < self.T_grid[0]:
    #         # print('inside less loop')
    #         it = np.searchsorted(self.T_grid, T, side='left')
    #         return self.Cp_grid[it - 1]*T
    #         # return self.Cp_grid[0]*T
    #     max_index = len(self.T_grid) - 1
    #     if T > self.T_grid[max_index]:
    #         # print('inside more loop')
    #         return self.H_grid[max_index] + self.Cp_grid[max_index] * (T - self.T_grid[max_index])
    #     else:
    #         it = np.searchsorted(self.T_grid, T, side='left')
    #         T0 = self.T_grid[it-1]
    #         T1 = self.T_grid[it]
    #         H_0 = self.H_grid[it-1]
    #         H_1 = self.H_grid[it]
    #         # print('inside general loop', it, T0, T1, H_0, H_1)
    #         return self.interpolated_value(T, T0, T1, H_0, H_1)

    # def get_mu(self):
    #     return self.mu
    #
    # def get_dH0(self):
    #     return self.dH0
    #
    #
    # def calculate_Cp_T_base(self):
    #     Cp_298 = self.Cp_From_Terra(self.T0)
    #     H_298 = self.H_From_Terra(self.T0)
    #     Cp_T_base = (H_298 - Cp_298*(self.T0-self.T_base)/2) / (self.T_base + (self.T0-self.T_base)/2)
    #     return Cp_T_base


# class Component_Dispersed(Component):
#     def __init__(self, name, T_base, mu, rho, dH0, T_grid, Cp_grid, H_grid):
#         super().__init__(name, T_base, mu, rho, dH0, T_grid, Cp_grid, H_grid)
#         self.columns_mix = ['T [K]', 'Cp [Дж/кг-К]', 'H [Дж/кг]']
#         self.data = pd.DataFrame(data=np.array([self.T_grid, self.Cp_grid, self.H_grid]).transpose(), index=self.T_grid, columns=self.columns_mix)
#         # self.data.drop(index=0, inplace=True)
#         print('mixture dispersed component ', name, ' is initialized')
#
#     def get_rho(self):
#         return self.rho
#
#
# class Component_Gas(Component):
#     def __init__(self, name, T_base, mu, rho, dH0, T_grid, Cp_grid, H_grid, viscosity_grid, heat_conductivity_grid, diffusivity_grid):
#         super().__init__(name, T_base, mu, rho, dH0, T_grid, Cp_grid, H_grid)
#         self.viscosity_grid = viscosity_grid
#         self.heat_conductivity_grid = heat_conductivity_grid
#         self.diffusivity_grid = diffusivity_grid
#         self.columns_mix = ['T [K]', 'Cp [Дж/кг-К]', 'H [Дж/кг]', 'Mu_visc [kg/m-s]', 'Lambda [W/m-K]', 'D [m^2/s]']
#         self.data = pd.DataFrame(data=np.array([self.T_grid, self.Cp_grid, self.H_grid, self.viscosity_grid, self.heat_conductivity_grid, self.diffusivity_grid]).transpose(), index=self.T_grid, columns=self.columns_mix)
#         # self.data.drop(index=0, inplace=True)
#         print('mixture gas component ', name, ' is initialized')
#
#
#     def get_viscosity(self, T):
#         return self.get_value_from_grid(T, self.viscosity_grid)
#
#     def get_heat_conductivity(self, T):
#         return self.get_value_from_grid(T, self.heat_conductivity_grid)
#
#     def get_diffusivity(self, T):
#         return self.get_value_from_grid(T, self.diffusivity_grid)
