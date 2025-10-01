import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import yaml
from matplotlib.patches import Polygon
from decimal import Decimal

from . import output as ou

# Константы:
R0 = 8.31




def initialize_terra_component(name, df_props, T0, T_base, dT_phase_transition, T_first, T_last, dT):
    print('inside initialize_terra_component', name)
    if name in df_props.index:
        _mu = df_props['mu'][name]
        _rho = df_props['rho'][name]
        _sigma = df_props['sigma'][name]
        _eps_dk = df_props['eps_dk'][name]
        _delta = df_props['delta'][name]
        _dH0 = df_props['dH0'][name]
        _T_range = df_props['T_range'][name]
        _f_ranges = df_props['f_ranges'][name]
        print(_mu)
        return Component(name, _mu, _rho, _sigma, _eps_dk, _delta, _dH0, _T_range, _f_ranges, T0, T_base, dT_phase_transition, T_first, T_last, dT)
    else:
        print('component with name ', name, ' not found in the database')
        return None


def get_prop(df, str_start, x, y, quant, str_end):
    # print('data from get props: ', df)
    prop = str_start
    prop_points = ''
    n_points = 0
    for ind, item in enumerate(df[x]):
        if not np.isnan(df[y][item]):
            n_points += 1
            prop_points += ''.join((' ', str(df[x][item]), ' ', str(Decimal(str(df[y][item])).quantize(Decimal(quant)))))
            print(prop_points)
    prop += str(n_points)
    prop += prop_points
    prop += str_end
    return prop

def show_properties_of_component(component, show_plots, T_base):
    is_gas = component.is_gas
    print(component.name, ' is gas: ', is_gas)
    name = component.name
    df = component.data

    T_PT = component.phase_transition_temperatures
    dT_PT = component.dT_phase_transition
    polygons = component.polygons
    vertical_lines = []
    for T in T_PT:
        # vertical_lines.append(T - dT_PT/2)
        # vertical_lines.append(T + dT_PT/2)
        vertical_lines.append(T)

    x_T_label = r'$T\ (K)$'
    y_Cp_label = r'$Cp\ (\frac{Дж}{кг \cdot К})$'
    y_H_label = r'$H\ (\frac{Дж}{кг})$'

    if show_plots:
        df.plot(y='Cp_from_Terra [Дж/кг-К]', use_index=True, style='-', grid=True, label='Cp Terra', xlabel=x_T_label, ylabel=y_Cp_label, xlim=0, ylim=0)
        plt.savefig(''.join(('pic_', name, '_01_Cp_Terra.jpeg')), dpi=400, bbox_inches='tight')
        plt.close()
        df.plot(y='Cp [Дж/кг-К]', use_index=True, style='o-', grid=True, label='Cp', xlabel=x_T_label, ylabel=y_Cp_label, xlim=0, ylim=0)
        plt.savefig(''.join(('pic_', name, '_02_Cp_Fixed_.jpeg')), dpi=400, bbox_inches='tight')
        plt.close()
        df.plot(y=['Cp [Дж/кг-К]', 'Cp_from_Terra [Дж/кг-К]'], use_index=True, style=['o-', '-'], grid=True, label=['Cp', 'Cp Terra'], xlabel=x_T_label, ylabel=y_Cp_label, xlim=0, ylim=0)
        for line in vertical_lines:
            plt.axvline(line, color='black', linestyle='--', label=''.join(('граница полинома (ФП),\nT = ', str(line), ' K')))
        plt.legend(loc='best', fontsize='small')
        ax = plt.gca()
        for polygon in polygons:
            for subpolygon in polygon:
                # print(subpolygon)
                p = Polygon(subpolygon, facecolor='C1')
                ax.add_patch(p)
        plt.savefig(''.join(('pic_', name, '_03_Cp_Compare.jpeg')), dpi=400, bbox_inches='tight')
        plt.close()

        df.plot(y='H_from_Terra [Дж/кг]', use_index=True, style='-', grid=True, label='H Terra', xlabel=x_T_label, ylabel=y_H_label, xlim=0, ylim=0)
        plt.savefig(''.join(('pic_', name, '_04_H_Terra.jpeg')), dpi=400, bbox_inches='tight')
        plt.close()
        df.plot(y='H [Дж/кг]', use_index=True, style='o-', grid=True, label='H', xlabel=x_T_label, ylabel=y_H_label, xlim=0, ylim=0)
        plt.savefig(''.join(('pic_', name, '_05_H_Integral.jpeg')), dpi=400, bbox_inches='tight')
        plt.close()
        df.plot(y=['H [Дж/кг]', 'H_from_Terra [Дж/кг]'], use_index=True, style=['o-', '-'], grid=True, label=['H', 'H Terra'], xlabel=x_T_label, ylabel=y_H_label, xlim=0, ylim=0)
        for line in vertical_lines:
            plt.axvline(line, color='black', linestyle='--', label=''.join(('граница полинома (ФП),\nT = ', str(line), ' K')))
        plt.legend(loc='best', fontsize='small')
        plt.savefig(''.join(('pic_', name, '_06_H_Compare.jpeg')), dpi=400, bbox_inches='tight')
        plt.close()

    # Make output files - yaml for QUBIC gas and dispersed phase / txt for fulent UDF dispersed phase / jou for fluent gas phase
    return_value = ou.output_props(name, is_gas, T_base, df)
    # return return_df_for_dict.index.tolist()
    return return_value

class Component:
    def __init__(self, name, mu, rho, sigma, eps_dk, delta, dH0, T_range, f_ranges, T0, T_base, dT_phase_transition, T_first, T_last, dT):
        self.name = name
        self.mu = mu
        self.rho = rho
        self.sigma = sigma
        self.eps_dk = eps_dk
        self.delta = delta
        if self.rho == 0 and self.sigma != 0 and self.eps_dk != 0:
            self.is_gas = True
            print('component ', self.name, ' is gas phase')
        elif self.rho != 0 and self.sigma == 0 and self.eps_dk == 0:
            self.is_gas = False
            print('component ', self.name, ' is dispersed phase')
        else:
            print('component ', self.name, ' is undefined phase (gas or dispersed?), check input properties: rho, sigma, eps_dk: ', self.rho, self.sigma, self.eps_dk)
        self.dH0 = dH0
        self.T_range = T_range
        self.f_ranges = f_ranges
        self.T0 = T0
        self.T_base = T_base
        self.dT_phase_transition = dT_phase_transition
        self.T_first = T_first
        self.T_last = T_last
        self.dT = dT
        self.phase_transition = self.Phase_Transition()
        self.phase_transition_temperatures = self.phase_transition[0]
        self.phase_transition_heats = self.phase_transition[1]

        self.vector_Cp = np.vectorize(self.Cp_From_Terra)

        self.T_grid = np.arange(self.T_first, self.T_last + self.dT, self.dT)
        self.update_T_grid()
        # print('T_grid is: ', self.T_grid)

        self.Cp_grid = self.vector_Cp(T=self.T_grid)
        self.Cp_phase_transition = self.calculate_Cp_phase_transition()[0]
        self.polygons = self.calculate_Cp_phase_transition()[1]
        print('Cp phase transition is: ', self.Cp_phase_transition)
        self.Cp_Tbase = self.calculate_Cp_T_base()

        self.vector_Cp_corrected = np.vectorize(self.Cp_corrected)
        self.Cp_corrected_grid = self.vector_Cp_corrected(T=self.T_grid)

        self.H_grid = np.zeros(shape=(len(self.T_grid)))
        self.update_H_grid()


        self.vector_H = np.vectorize(self.H_From_Terra)
        self.H_grid_from_Terra = self.vector_H(T=self.T_grid)

        print('T_grid is: ')
        print(self.T_grid)

        if not self.is_gas:
            self.columns_mix = ['T [K]', 'Cp [Дж/кг-К]', 'H [Дж/кг]', 'Cp_from_Terra [Дж/кг-К]', 'H_from_Terra [Дж/кг]']
            self.data = pd.DataFrame(data=np.array([self.T_grid, self.Cp_corrected_grid, self.H_grid, self.Cp_grid, self.H_grid_from_Terra]).transpose(), index=self.T_grid, columns=self.columns_mix)
            print('dispersed component ', name, ' is initialized')
        else:
            self.columns_mix = ['T [K]', 'Cp [Дж/кг-К]', 'H [Дж/кг]', 'Cp_from_Terra [Дж/кг-К]', 'H_from_Terra [Дж/кг]', 'Mu_visc [kg/m-s]', 'Lambda [W/m-K]', 'D [m^2/s]']
            self.vector_viscosity = np.vectorize(self._viscosity_kinetic_theory)
            self.viscosity_grid = self.vector_viscosity(T=self.T_grid)
            self.vector_heat_conductivity = np.vectorize(self._heat_conductivity_kinetic_theory)
            self.heat_conductivity_grid = self.vector_heat_conductivity(T=self.T_grid)
            self.vector_diffusivity = np.vectorize(self._diffusivity_kinetic_theory)
            self.diffusivity_grid = self.vector_diffusivity(T=self.T_grid)
            self.data = pd.DataFrame(data=np.array([self.T_grid, self.Cp_corrected_grid, self.H_grid, self.Cp_grid, self.H_grid_from_Terra, self.viscosity_grid, self.heat_conductivity_grid, self.diffusivity_grid]).transpose(), index=self.T_grid, columns=self.columns_mix)
            print('gas component ', name, ' is initialized')

        self.data.drop(index=0, inplace=True)
        self.writer = pd.ExcelWriter(''.join(('out_', self.name, '.xlsx')), engine="xlsxwriter")
        self.data.to_excel(self.writer, index=False, sheet_name='props')
        self.writer.close()


    def update_T_grid(self):
        T_append = np.array([self.T_range[0], self.T_range[-1]])
        ind_exclude_because_of_min_max_limits = np.nonzero((T_append <= self.T_first) | (T_append >= self.T_last))[0]
        T_append = np.delete(T_append, ind_exclude_because_of_min_max_limits)

        indexes_between_to_remove = []
        for T in self.phase_transition_temperatures:
            if (T > (self.T_first + self.dT_phase_transition / 2)) and (T < (self.T_last - self.dT_phase_transition / 2)):
                T_before = T - self.dT_phase_transition / 2
                T_after = T + self.dT_phase_transition / 2
                T_append = np.append(T_append, np.array([T_before, T_after]))
                indices_between = (np.nonzero((self.T_grid >= T_before) & (self.T_grid <= T_after))[0]).tolist()
                indexes_between_to_remove += indices_between
        self.T_grid = np.delete(self.T_grid, indexes_between_to_remove)
        self.T_grid = np.append(self.T_grid, T_append)
        self.T_grid = np.unique(self.T_grid)
        self.T_grid.sort()
        indexes_check = np.where(np.in1d(self.T_grid, T_append))[0]
        indexes_remove = []
        for i in indexes_check:
            if i != 0 and i != len(self.T_grid - 1):
                T_before = self.T_grid[i - 1]
                T_after = self.T_grid[i + 1]
                T_current = self.T_grid[i]
                if abs(T_current - T_before) < self.dT / 2:
                    indexes_remove.append(i - 1)
                if abs(T_after - T_current) < self.dT / 2:
                    indexes_remove.append(i + 1)
        self.T_grid = np.delete(self.T_grid, indexes_remove)
        return None












    def S(self, x, f1, f2, f3, f5, f6, f7):
        return (f1 + f2 + f2 * math.log(x) - f3 * pow(x, -2) + 2 * f5 * x + 3 * f6 * pow(x, 2) + 4 * f7 * pow(x, 3)) / self.mu * 1000

    def H(self, x, f2, f3, f4, f5, f6, f7):
        return (f2*x - 2*f3*pow(x, -1) - f4 + f5*pow(x,2) + 2*f6*pow(x, 3) + 3*f7*pow(x, 4))*10000 / self.mu * 1000

    def Cp(self, x, f2, f3, f5, f6, f7):
        return (f2 + 2*f3*pow(x, -2) + 2*f5*x + 6*f6*pow(x, 2) + 12*f7*pow(x, 3)) / self.mu * 1000


    def Phase_Transition(self):
        T_pt = self.T_range[1:-1]
        Q_pt = np.zeros(shape=(len(T_pt)))
        # print('inside')
        for i, Q in enumerate(Q_pt):
            f_before = self.f_ranges[i]
            f_after = self.f_ranges[i + 1]
            x = T_pt[i] * 0.0001
            S_before = self.S(x, f_before[0], f_before[1], f_before[2], f_before[4], f_before[5], f_before[6])
            S_after = self.S(x, f_after[0], f_after[1], f_after[2], f_after[4], f_after[5], f_after[6])
            dS = S_after - S_before
            Q_pt[i] = dS * T_pt[i]
            # print('Phase transition point: ', T_pt[i], Q_pt[i])
        return T_pt, Q_pt

    def Cp_From_Terra(self, T):
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
        # print(i)
        x = 0.0001*T
        return self.Cp(x, self.f_ranges[i, 1], self.f_ranges[i, 2], self.f_ranges[i, 4], self.f_ranges[i, 5], self.f_ranges[i, 6])

    def H_From_Terra(self, T):
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
        # print(i)
        x = 0.0001*T
        return self.H(x, self.f_ranges[i, 1], self.f_ranges[i, 2], self.f_ranges[i, 3], self.f_ranges[i, 4], self.f_ranges[i, 5], self.f_ranges[i, 6])

    def calculate_Cp_T_base(self):
        Cp_298 = self.Cp_From_Terra(self.T0)
        H_298 = self.H_From_Terra(self.T0)
        Cp_T_base = (H_298 - Cp_298*(self.T0-self.T_base)/2) / (self.T_base + (self.T0-self.T_base)/2)
        return Cp_T_base

    def calculate_Cp_phase_transition(self):
        Cp_PT = np.zeros(shape=(len(self.phase_transition_temperatures)))
        polygons = []
        for i, T in enumerate(self.phase_transition_temperatures):
            T_b = T - self.dT_phase_transition / 2  # температура до
            T_a = T + self.dT_phase_transition / 2  # температура после
            ind_before = np.nonzero(self.T_grid == T_b)[0][0]
            ind_after = np.nonzero(self.T_grid == T_a)[0][0]
            T_bb = self.T_grid[ind_before - 1]  # температура до до
            T_aa = self.T_grid[ind_after + 1]  # температура после после
            Cp_b = self.Cp_grid[ind_before]
            Cp_a = self.Cp_grid[ind_after]
            # only for polygons - Cp_bb Cp_aa
            Cp_bb = self.Cp_grid[ind_before - 1]
            Cp_aa = self.Cp_grid[ind_after + 1]
            dT_bb = T_b - T_bb
            dT_ba = T_a - T_b
            dT_aa = T_aa - T_a
            Cp_PT[i] = (Cp_b * dT_bb / 2 + (Cp_b + Cp_a) * dT_ba / 2 + Cp_a * dT_aa / 2 + self.phase_transition_heats[i]) / (dT_bb / 2 + dT_ba + dT_aa / 2)
            # print('calculated ind T Cp: ', ind_before, ind_after, T_bb, T_b, T_a, T_aa, Cp_b, Cp_a, Cp_PT)
            # Now calculate reactangle coodinates to show on plots
            polygons.append([np.array([[T_bb, Cp_bb], [T_b, Cp_PT[i]], [T_b, Cp_b], [T_bb, Cp_bb]]),
                             np.array([[T_b, Cp_b], [T_b, Cp_PT[i]], [T_a, Cp_PT[i]], [T_a, Cp_a], [T_b, Cp_b]]),
                             np.array([[T_a, Cp_a], [T_a, Cp_PT[i]], [T_aa, Cp_aa], [T_a, Cp_a]])])
        return Cp_PT, polygons

    def interpolated_value(self, x, x0, x1, y0, y1):
        return y0 + (y1 - y0) / (x1 - x0) * (x - x0)

    # Cp with all corrections to account for phase transition and to make enthalpy from this Cp integration exactly equal terra enthalpy values
    def Cp_corrected(self, T):
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
            if T <= self.T_base:
                return self.Cp_Tbase
            if T > self.T_base:
                return self.interpolated_value(T, self.T_base, self.T0, self.Cp_Tbase, self.Cp_From_Terra(T))
            # T = T_range[0]
        if (T >= (self.T_range[i + 1] - self.dT_phase_transition / 2)) and (i <= (len(self.T_range) - 3)):
            # print('first', i, T)
            return self.Cp_phase_transition[i]
        elif (T <= (self.T_range[i] + self.dT_phase_transition / 2)) and (i != 0):
            # print('second', i, T)
            return self.Cp_phase_transition[i - 1]
        else:
            # print('pass', i, T)
            return self.Cp_From_Terra(T)

    def update_H_grid(self):
        H_prev = 0
        for i, value in enumerate(self.H_grid):
            if i == 0:
                dH = self.Cp_corrected_grid[i]*self.T_grid[i]
            else:
                Cp_prev = self.Cp_corrected_grid[i-1]
                Cp_cur = self.Cp_corrected_grid[i]
                Cp_mean = (Cp_prev + Cp_cur)/2
                dT = self.T_grid[i] - self.T_grid[i-1]
                dH = Cp_mean*dT
            H_current = H_prev + dH
            self.H_grid[i] = H_current
            H_prev = H_current
        return None

    # KINETIC THEORY PARAMETERS CALCULATION:
    # reference temperature for diffusivity
    def _T_D(self, T):
        if T < self.T_base:
            T = self.T_base
        # print('theta T: ', T)
        return T / self.eps_dk

    # collision integral for diffusivity
    def _theta_diffusivity(self, T):
        return 1 / pow(self._T_D(T), 0.145) + 1 / pow((self._T_D(T) + 0.5),2)

    # collision integral for viscosity
    def _theta_viscosity(self, T):
        return 1.14 * self._theta_diffusivity(T)

    # viscocity [kg/m-s]
    def _viscosity_kinetic_theory(self, T):
        return 2.67e-6 * pow(self.mu * T, 0.5) / (pow(self.sigma, 2) * self._theta_viscosity(T))

    # heat conductivity [W/m-K]
    def _heat_conductivity_kinetic_theory(self, T):
        return 15/4 * R0 / (self.mu * 0.001) * self._viscosity_kinetic_theory(T) * (4/15 * self.get_Cp(T)*self.mu*0.001/R0 + 1/3)

    # self-diffusivity [m^2/s]
    def _diffusivity_kinetic_theory(self, T):
        return 0.00186 * pow(pow(T, 3) * (1/self.mu + 1/self.mu), 0.5) / (pow(self.sigma, 2) * self._theta_diffusivity(T))


    def get_Cp(self, T):
        if T < min(self.T_grid):
            return self.Cp_corrected_grid[0]
        if T > max(self.T_grid):
            return self.Cp_corrected_grid[len(self.Cp_corrected_grid)-1]
        else:
            it = np.searchsorted(self.T_grid, T, side='left')
            T0 = self.T_grid[it-1]
            T1 = self.T_grid[it]
            Cp_0 = self.Cp_corrected_grid[it-1]
            Cp_1 = self.Cp_corrected_grid[it]
            return self.interpolated_value(T, T0, T1, Cp_0, Cp_1)

    def get_H(self, T):
        if T < self.T_grid[0]:
            return self.Cp_corrected_grid[0]*T
        max_index = len(self.T_grid) - 1
        if T > self.T_grid[max_index]:
            return self.H_grid[max_index] + self.Cp_corrected_grid[max_index] * (T - self.T_grid[max_index])
        else:
            it = np.searchsorted(self.T_grid, T, side='left')
            T0 = self.T_grid[it-1]
            T1 = self.T_grid[it]
            H_0 = self.H_grid[it-1]
            H_1 = self.H_grid[it]
            return self.interpolated_value(T, T0, T1, H_0, H_1)

    def get_viscosity(self, T):
        if T < self.T_base:
            return self.viscosity_grid[1]
        if T > max(self.T_grid):
            return self.viscosity_grid[len(self.viscosity_grid)-1]
        else:
            it = np.searchsorted(self.T_grid, T, side='left')
            T0 = self.T_grid[it-1]
            T1 = self.T_grid[it]
            y_0 = self.viscosity_grid[it-1]
            y_1 = self.viscosity_grid[it]
            return self.interpolated_value(T, T0, T1, y_0, y_1)

    def get_heat_conductivity(self, T):
        if T < self.T_base:
            return self.heat_conductivity_grid[1]
        if T > max(self.T_grid):
            return self.heat_conductivity_grid[len(self.heat_conductivity_grid)-1]
        else:
            it = np.searchsorted(self.T_grid, T, side='left')
            T0 = self.T_grid[it-1]
            T1 = self.T_grid[it]
            y_0 = self.heat_conductivity_grid[it-1]
            y_1 = self.heat_conductivity_grid[it]
            return self.interpolated_value(T, T0, T1, y_0, y_1)

    def get_diffusivity(self, T):
        if T < self.T_base:
            return self.diffusivity_grid[1]
        if T > max(self.T_grid):
            return self.diffusivity_grid[len(self.diffusivity_grid)-1]
        else:
            it = np.searchsorted(self.T_grid, T, side='left')
            T0 = self.T_grid[it-1]
            T1 = self.T_grid[it]
            y_0 = self.diffusivity_grid[it-1]
            y_1 = self.diffusivity_grid[it]
            return self.interpolated_value(T, T0, T1, y_0, y_1)

    def get_mu(self):
        return self.mu

    def get_rho(self):
        return self.rho

    def get_dH0(self):
        return self.dH0

    def get_name(self):
        return self.name






