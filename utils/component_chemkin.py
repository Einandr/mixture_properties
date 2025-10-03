import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mendeleev import element
import re
import math
import os

# Константы:
R0 = 8.314       # [Дж / моль-К] универсальная газовая постоянная
cal = 4184       # [кДж] термическая калория, используемая в Ansys Chemkin


class Component:
    def __init__(self, name, date, formula, phase, T_low, T_mid, T_high, atomic, a1_low, a2_low, a3_low, a4_low, a5_low, a6_low, a7_low, a1_high, a2_high, a3_high, a4_high, a5_high, a6_high, a7_high, dT, T_base, T0, eps_dk, sigma, T_last=None):
        self.name = name
        self.date = date
        self.formula = formula
        self.phase = phase
        self.T_low = T_low
        self.T_mid = T_mid
        self.T_high = T_high
        self.T_last = T_last
        self.atomic = atomic
        self.T_range = np.array([T_low, T_mid, T_high])
        self.a = np.array([[a1_low, a2_low, a3_low, a4_low, a5_low, a6_low, a7_low], [a1_high, a2_high, a3_high, a4_high, a5_high, a6_high, a7_high]])
        self.dT = dT                                                        # интервал с которым будем создавать сетку значений по температуре
        self.T_base = T_base                                                # [K] температура для вычисления крайнего нижнего значения Ср с целью сведения в единое значение полной энтальпии при расчете через коэффициенты и через интеграл
        self.T0 = T0                                                        # должна быть 298.15 K - температура для вычисления энтальпии образования

        print('check formula:', formula)
        self.sp_list = (re.sub(r'(\d+[.,]?)', r'\1 ', self.formula))        # выделяем границу молей предыдущего компонента и имени следующего
        self.sp_list = (re.sub(r'([a-zA-Z]+)', r'\1 ', self.sp_list))       # выделяем границу компонента и его количества молей
        self.sp_list = re.split(r'\s', self.sp_list)                        # разбиваем строку по пробелам
        self.sp_list = [x for x in self.sp_list if x]

        self.sp_components = [x for ind, x in enumerate(self.sp_list) if ((ind % 2) != 1)]
        self.sp_moles = [x for ind, x in enumerate(self.sp_list) if ((ind % 2) != 0)]
        self.mu = 0                                                         # [г/моль]
        print('test sp_components:', self.sp_components, self.sp_moles)
        for ind, x in enumerate(self.sp_components):
            self.el = element(x)
            self.moles = float(self.sp_moles[ind])
            self.mu += self.el.atomic_weight * self.moles
        self.R = R0 / self.mu * 1000                                        # [Дж/(кг-К)]
        self.path_chemkin = 'chemkin'

        self.dH0 = self.H0_298_J()
        # self.T_grid = np.arange(self.T_low, self.T_high + self.dT, self.dT)
        self.T_grid = self.create_T_grid()
        self.Cp_Tbase = self.calculate_Cp_T_base()
        print('T_low, T_mid T_high, dT', self.T_low, self.T_mid, self.T_high, self.dT)
        print('self T_grid', self.T_grid)
        self.update_T_grid()
        print('T_bsase, T_mid:', self.T_base, self.T_mid)
        print('self T_grid', self.T_grid)

        self.vector_Cp = np.vectorize(self.Cp, excluded=['units'])
        self.Cp_grid = self.vector_Cp(T=self.T_grid, units='J_kg_K')
        print('self Cp_grid', self.Cp_grid)

        self.H_grid = np.zeros(shape=(len(self.T_grid)))
        self.update_H_grid()
        print('self H_grid', self.H_grid)

        self.rho = 0    # Заглушка для газа
        # Значения по умолчанию для воздуха для транспортных свойств
        self.eps_dk = eps_dk
        self.sigma = sigma

        self.columns_mix = ['T [K]', 'Cp [Дж/кг-К]', 'H [Дж/кг]', 'Mu_visc [kg/m-s]', 'Lambda [W/m-K]', 'D [m^2/s]']
        self.vector_viscosity = np.vectorize(self._viscosity_kinetic_theory)
        self.viscosity_grid = self.vector_viscosity(T=self.T_grid)
        self.vector_heat_conductivity = np.vectorize(self._heat_conductivity_kinetic_theory)
        self.heat_conductivity_grid = self.vector_heat_conductivity(T=self.T_grid)
        self.vector_diffusivity = np.vectorize(self._diffusivity_kinetic_theory)
        self.diffusivity_grid = self.vector_diffusivity(T=self.T_grid)
        self.data = pd.DataFrame(data=np.array([self.T_grid, self.Cp_grid, self.H_grid, self.viscosity_grid, self.heat_conductivity_grid, self.diffusivity_grid]).transpose(), index=self.T_grid, columns=self.columns_mix)
        print('gas component ', name, ' is initialized')
        print(self.data)

        os.makedirs(self.path_chemkin, exist_ok=True)
        os.chdir(self.path_chemkin)
        # self.data.drop(index=0, inplace=True)
        self.writer = pd.ExcelWriter(''.join(('out_', self.name, '.xlsx')), engine="xlsxwriter")
        self.data.to_excel(self.writer, index=False, sheet_name='props')
        self.writer.close()
        os.chdir('..')

    def create_T_grid(self):
        rounded_T_low = math.ceil(self.T_low / 100) * 100       # Округляем T_low до ближайшего целого, кратного 100, но не меньше T_low
        if rounded_T_low < self.T_low:
            rounded_T_low += 100
        T_upper_bound = self.T_high                                 # Определяем верхнюю границу температурной сетки
        if hasattr(self, 'T_last') and self.T_last is not None:
            T_upper_bound = self.T_last
        T_grid = np.arange(rounded_T_low, T_upper_bound + self.dT, self.dT)   # Создаем массив температур, начиная с rounded_T_low и заканчивая T_high
        print('before deletion T_low:', T_grid)
        T_grid = T_grid[1:]                                     # Удаляем первое значение - это либо T_low, либо rounded_T_low
        T_grid = np.insert(T_grid, 0, self.T_low)               # Добавляем обратно T_low
        T_grid = T_grid[T_grid <= T_upper_bound]                # Убедимся, что значения не превышают T_upper_bound
        if T_grid[-1] != T_upper_bound:                         # Проверяем, содержится ли T_upper_bound в массиве T_grid
            T_grid = np.append(T_grid, T_upper_bound)
        print('after checking T_grid:', T_grid)
        return T_grid

    def interpolated_value(self, x, x0, x1, y0, y1):
        # print('interpolating', x, x0, x1, y0, y1)
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

    def calculate_Cp_T_base(self):
        Cp_298 = self.Cp(self.T0, 'J_kg_K')
        H_298 = self.H0_298_J()
        Cp_T_base = (H_298 - Cp_298*(self.T0-self.T_base)/2) / (self.T_base + (self.T0-self.T_base)/2)
        return Cp_T_base

    def Cp(self, T, units):         # Эта функция используется для заполнения массива по запрашиваемой температуре
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

        elif T <= self.T_range[0]:
            if T <= self.T_base:
                return self.Cp_Tbase
            if T > self.T_base:
                if units == 'kcal_mole':
                    Cp_T_range0 =  self._Cp_kcal_mole(T, self.a[i, 0], self.a[i, 1], self.a[i, 2], self.a[i, 3], self.a[i, 4])
                elif units == 'J_kg_K':
                    Cp_T_range0 = self._Cp_J_kg_K(T, self.a[i, 0], self.a[i, 1], self.a[i, 2], self.a[i, 3], self.a[i, 4])
                else:
                    raise ValueError('неверные единицы измерения Ср')
                return self.interpolated_value(T, self.T_base, self.T0, self.Cp_Tbase, Cp_T_range0)
        if units == 'kcal_mole':
            return self._Cp_kcal_mole(T, self.a[i, 0], self.a[i, 1], self.a[i, 2], self.a[i, 3], self.a[i, 4])
        elif units == 'J_kg_K':
            return self._Cp_J_kg_K(T, self.a[i, 0], self.a[i, 1], self.a[i, 2], self.a[i, 3], self.a[i, 4])
        else:
            raise ValueError('неверные единицы измерения Ср')

    def update_T_grid(self):
        self.T_grid = np.append(self.T_grid, [self.T_mid, self.T_base])
        # добавляем T_high на случай если T_last > T_high и мы его проскочили при создании T_grid
        if hasattr(self, 'T_last') and self.T_last is not None and self.T_last > self.T_high:
            self.T_grid = np.append(self.T_grid, self.T_high)
        self.T_grid = np.unique(self.T_grid)
        self.T_grid.sort()
        return None

    def update_H_grid(self):
        H_prev = 0
        for i, value in enumerate(self.H_grid):
            if i == 0:
                dH = self.Cp_grid[i]*self.T_grid[i]
            else:
                Cp_prev = self.Cp_grid[i-1]
                Cp_cur = self.Cp_grid[i]
                Cp_mean = (Cp_prev + Cp_cur)/2
                dT = self.T_grid[i] - self.T_grid[i-1]
                dH = Cp_mean*dT
            H_current = H_prev + dH
            self.H_grid[i] = H_current
            H_prev = H_current
        return None


    # KINETIC THEORY PARAMETERS CALCULATION: (ДУБЛИРОВАНИЕ COMPONENT_TERRA)

    # reference temperature for diffusivity
    def _T_D(self, T):
        if T < self.T_base:
            T = self.T_base
        # print('theta T: ', T)
        return T / self.eps_dk

    # collision integral for diffusivity
    def _theta_diffusivity(self, T):
        return 1 / pow(self._T_D(T), 0.145) + 1 / pow((self._T_D(T) + 0.5), 2)

    # collision integral for viscosity
    def _theta_viscosity(self, T):
        return 1.14 * self._theta_diffusivity(T)

    # viscocity [kg/m-s]
    def _viscosity_kinetic_theory(self, T):
        return 2.67e-6 * pow(self.mu * T, 0.5) / (pow(self.sigma, 2) * self._theta_viscosity(T))

    # heat conductivity [W/m-K]
    def _heat_conductivity_kinetic_theory(self, T):
        return 15 / 4 * R0 / (self.mu * 0.001) * self._viscosity_kinetic_theory(T) * (4 / 15 * self.get_Cp(T) * self.mu * 0.001 / R0 + 1 / 3)

    # self-diffusivity [m^2/s]
    def _diffusivity_kinetic_theory(self, T):
        return 0.00186 * pow(pow(T, 3) * (1 / self.mu + 1 / self.mu), 0.5) / (pow(self.sigma, 2) * self._theta_diffusivity(T))


    # Эти методы должны быть с теми же названиями, что для component_terra
    def get_Cp(self, T):
        return self.get_value_from_grid(T, self.Cp_grid)

    def get_H(self, T):
        if T < self.T_grid[0]:
            return self.Cp_grid[0]*T
        max_index = len(self.T_grid) - 1
        if T > self.T_grid[max_index]:
            return self.H_grid[max_index] + self.Cp_grid[max_index] * (T - self.T_grid[max_index])
        else:
            it = np.searchsorted(self.T_grid, T, side='left')
            T0 = self.T_grid[it-1]
            T1 = self.T_grid[it]
            H_0 = self.H_grid[it-1]
            H_1 = self.H_grid[it]
            return self.interpolated_value(T, T0, T1, H_0, H_1)

    def get_mu(self):
        return self.mu

    def get_dH0(self):
        return self.dH0

    def get_rho(self):
        return self.rho

    def get_viscosity(self, T):
        return self.get_value_from_grid(T, self.viscosity_grid)

    def get_heat_conductivity(self, T):
        return self.get_value_from_grid(T, self.heat_conductivity_grid)

    def get_diffusivity(self, T):
        return self.get_value_from_grid(T, self.diffusivity_grid)


