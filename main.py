import math
import dearpygui.dearpygui as dpg

MIN_WIDTH = 1000
MIN_HEIGHT = 700

R_inner = 0.305 / 2      # Внутренний радиус (м)
R_outer = 2 * R_inner # Внешний радиус (м)
wall_thickness = R_outer - R_inner

T_initial_C = 15.0   # Начальная температура стенки (°C)
T_air_C = 20.0       # Температура окружающего воздуха (°C)

# Перевод температур в Кельвины
T_initial_K = T_initial_C + 273.15
T_air_K = T_air_C + 273.15
T_gas_K = 1075.0     # Температура пороховых газов (K) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Свойства материала (сталь)
DENSITY = 7850       # Плотность (кг/м^3)
LAMBDA = 34          # Теплопроводность (Вт/(м*К))
C_HEAT = 496         # Удельная теплоемкость (Дж/(кг*К)

# Коэффициенты теплоотдачи
ALPHA_GAS = 19600    # Коэфф. теплоотдачи с пороховым газом (Вт/(м^2*К)) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ALPHA_AIR = 9        # Коэфф. теплоотдачи с окружающим воздухом (Вт/(м^2*К))

SHOT_DURATION = 0.020 # Время выстрела (с)
K_TIME_STEPS = 20     # Количество шагов по времени

N_NODES = 51  # Количество узлов по радиусу (50 шагов)
H_MIN = 1e-6  # Минимальный шаг на внутренней поверхности (м)
num_steps_h = N_NODES - 1

b_cooling = 0.0037 # Параметр экспоненциального охлаждения газов

# Серия выстрелов
N_SHOTS = 1
FIRE_RATE_PER_MIN = 2.0  # 2 выстрела в минуту
INTERVAL_BETWEEN_SHOTS = 60.0 / FIRE_RATE_PER_MIN  # сек (включая сам выстрел)
COOL_DOWN_INTERVAL = INTERVAL_BETWEEN_SHOTS - SHOT_DURATION

def geometric_progression_ratio(L, n, h_min, eps=1e-6):
    """
    Определяет коэффициент q геометрической прогрессии.
    Сумма слоев стремится к L.

    Аргументы:
        L (float): общая толщина стенки (м)
        n (int): количество слоев
        h_min (float): минимальный шаг на внутренней поверхности (м)

    Возвращает:
        float: коэффициент q
    """
    if not (L > 0 and n > 0 and h_min > 0 and h_min < L):
        raise ValueError("Неверные входные параметры.")

    q_low, q_high = 1.0, 2.0
    while h_min * (q_high**n - 1) / (q_high - 1) < L:
        q_high *= 2

    while abs(q_high - q_low) > eps:
        q = 0.5 * (q_low + q_high)
        total_thickness = h_min * (q**n - 1) / (q - 1)
        if total_thickness > L:
            q_high = q
        else:
            q_low = q

    return q

def solve_barrel_heating(total_time, time_steps, h, r, initial_temp, gas_temp, air_temp, lambda_steel, c_steel, alpha_in, alpha_out, rho, b_cooling, flag):
    """
    Рассчитывает нагрев цилиндрической стенки с помощью неявного метода конечных разностей.

    Аргументы:
        total_time  (float): Общее время процесса в секундах.
        time_steps (int): Количество шагов по времени.
        h  (list): Cписок шагов по толщине (м).
        r  (list): Cписок радиусов (м).
        initial_temp  (list): Начальное распределение температуры по узлам (K).
        gas_temp  (float): Температура порозовых газов (K).
        air_temp  (float): Температура внешнего окружающего воздуха (K).
        lambda_steel  (float): Теплопроводность стали (Вт/(м*К)).
        c_steel  (float): Удельная теплоемкость стали (Дж/(кг*К)).
        alpha_in  (float): Коэффициент теплоотдачи на внутренней поверхности (Вт/(м^2*К)).
        alpha_out  (float): Коэффициент теплоотдачи на внешней поверхности (Вт/(м^2*К)).
        rho  (float): Плотность стали (кг/м^3).

    Возвращает:
        list: Список списков, содержащий температуру в каждом узле для каждого шага по времени.
    """
    N = len(r) #Количество узлов по толщине
    n = N - 1 #Число шагов по толщине
    delta_t = total_time / time_steps #Шаг по времени
    #Узлы по времени
    t = [0 + i*delta_t for i in range(time_steps + 1)]

    # Инициализация матрицы температур начальным условием
    U = [list(initial_temp)]
    
    def Lambda(a):
        """ Функция коэффициента теплопроводности материала трубы Вт/(м*°К) """
        return lambda_steel
    
    def CFunc(a):
        """ Функция - теплоемкость от температуры Дж/(кг*К)  """
        return c_steel
    
    def AlphaIn(t):
        """ Функция коэффициента теплоотдачи газа от времени на внутренней поверхности """
        if flag:
            return alpha_out + (alpha_in - alpha_out) * math.exp(-t / b_cooling*0.8)
        else:
            return alpha_in
    
    def Ug(t):
        """ Функция внутренней температуры газа от времени """
        if flag:
            # return air_temp + (initial_temp[0] - air_temp) * math.exp(-t / b_cooling*0.8)
            return air_temp + (gas_temp - air_temp) * math.exp(-t / b_cooling*0.8)
        else:
            return gas_temp

    #Расчёт
    for i in range(1, time_steps):
        #Создаем пустую строку температур
        U.append([None for i in range(N)])
        #Создаем строки прогоночных коэффициентов
        Alpha = [None for j in range(n)]
        Beta = [None for j in range(n)]
        #Прогоночные коэффициенты для первого узла
        Alpha[0] = Lambda(U[i-1][0])/(Lambda(U[i-1][0]) + h[0]*AlphaIn(t[i]))
        Beta[0] = h[0]*AlphaIn(t[i])*Ug(t[i])/(Lambda(U[i-1][0]) + h[0]*AlphaIn(t[i]))
        #Прогоночные коэффициенты для других точек
        for j in range(1, n):
            AJ = (Lambda(U[i - 1][j + 1]) + Lambda(U[i - 1][j]))/h[j]/(h[j] + h[j - 1]) + Lambda(U[i - 1][j])/r[j]/(h[j] + h[j - 1])
            CJ = (Lambda(U[i - 1][j]) + Lambda(U[i - 1][j - 1]))/h[j - 1]/(h[j] + h[j - 1]) - Lambda(U[i - 1][j])/r[j]/(h[j] + h[j - 1])
            BJ = CFunc(U[i - 1][j])*rho/delta_t + (Lambda(U[i - 1][j + 1]) + Lambda(U[i - 1][j]))/h[j]/(h[j] + h[j - 1]) + (Lambda(U[i - 1][j]) + Lambda(U[i - 1][j - 1]))/h[j - 1]/(h[j] + h[j - 1])
            FJ = -CFunc(U[i - 1][j])*rho*U[i - 1][j]/delta_t
            Alpha[j] = AJ/(BJ - CJ*Alpha[j - 1])
            Beta[j] = (CJ*Beta[j - 1] - FJ)/(BJ - CJ*Alpha[j - 1])
        #Температура на наружной поверхности цилиндра
        U[i][n] = (Lambda(U[i - 1][n])*Beta[n - 1] + h[n - 1]*alpha_out*air_temp)/(h[n - 1]*alpha_out + Lambda(U[i - 1][n])*(1 - Alpha[n - 1]))
        #Температура в других точках
        for j in range(n - 1, -1, -1):
            U[i][j] = Alpha[j]*U[i][j + 1] + Beta[j]
    #Возвращаем значение
    return U

def calc(T_gass, T_interes, alphaa, time_interes):
    strog = 2
    t_strog = 0.05
    T_interes_K = T_interes + 273.15
    dict1 = {}
    
    # Расчет коэффициента для геометрической прогрессии размеров шагов
    q_ratio = geometric_progression_ratio(wall_thickness, num_steps_h, H_MIN)
    
    # Создание списка размеров шагов и радиусов
    h_radial_steps = [H_MIN * q_ratio**i for i in range(num_steps_h)]
    radii_nodes = [R_inner]
    for i in range(num_steps_h):
        radii_nodes.append(radii_nodes[-1] + h_radial_steps[i])
    # radii_nodes = np.linspace(R_inner, R_outer, N_NODES)
    
    cur_temperatures = [T_initial_K] * N_NODES
    max_temp = 0
    
    rate_min = FIRE_RATE_PER_MIN
    rate_max = FIRE_RATE_PER_MIN * 5
    last_rate = 0
    solution_found = False  # флаг найденного решения
    
    max_iterations = 50  # защита от бесконечного цикла
    iteration_count = 0
    t_time = 0
    
    while iteration_count < max_iterations:
        iteration_count += 1
        
        temp_diff = abs(T_interes_K - max_temp)
        time_diff = abs(t_time - time_interes)
        
        if temp_diff <= strog and time_diff <= time_interes*t_strog:
            solution_found = True
            break
            
        if (rate_max - rate_min) < 0.01:
            break
        
        all_T = []
        ttime = [0]
        t_time = 0
        rate = (rate_max + rate_min) / 2
        interval_between_shots = 60.0 / rate
        cool_down_interval = interval_between_shots - SHOT_DURATION
        cur_temperatures = [T_initial_K] * N_NODES
        temperature_exceeded = False
        cur_temperatures_to_return = []
        ct = 0
        
        while t_time + interval_between_shots <= time_interes and not temperature_exceeded:
            ct += 1
            # Выстрел
            U_shot = solve_barrel_heating(
                total_time=SHOT_DURATION,
                time_steps=K_TIME_STEPS,
                h=h_radial_steps,
                r=radii_nodes,
                initial_temp=cur_temperatures,
                gas_temp=T_gass,
                air_temp=T_air_K,
                lambda_steel=LAMBDA,
                c_steel=C_HEAT,
                alpha_in=alphaa,
                alpha_out=ALPHA_AIR,
                rho=DENSITY,
                b_cooling=b_cooling,
                flag=False,
            )
            
            time_step_shot = SHOT_DURATION / K_TIME_STEPS
            for i, temp_step in enumerate(U_shot):
                current_time = t_time + (i + 1) * time_step_shot
                if temp_step[0] > T_interes_K + strog:
                    temperature_exceeded = True
                all_T.append(temp_step[0])
                ttime.append(current_time)
            
            t_time += SHOT_DURATION
            cur_temperatures = U_shot[-1]
            max_temp = cur_temperatures[0]
            cur_temperatures_to_return = cur_temperatures
        
            # Охлаждение между выстрелами
            U_cool = solve_barrel_heating(
                total_time=cool_down_interval,
                time_steps=K_TIME_STEPS,
                h=h_radial_steps,
                r=radii_nodes,
                initial_temp=cur_temperatures,
                gas_temp=T_gass,
                air_temp=T_air_K,
                lambda_steel=LAMBDA,
                c_steel=C_HEAT,
                alpha_in=alphaa,
                alpha_out=ALPHA_AIR,
                rho=DENSITY,
                b_cooling=b_cooling,
                flag=True,
            )
            
            time_step_cool = cool_down_interval / K_TIME_STEPS
            for i, temp_step in enumerate(U_cool):
                current_time = t_time + (i + 1) * time_step_cool
                if temp_step[0] > T_interes_K + strog:
                    temperature_exceeded = True
                all_T.append(temp_step[0])
                ttime.append(current_time)
            
            t_time += cool_down_interval
            cur_temperatures = U_cool[-1]
            max_temp = max([max(i) for i in U_cool])
            for i in U_cool:
                if (i[0]) == max_temp:
                    cur_temperatures_to_return = i
        
        print(f"Скорострельность: {rate_min:.2f} - {rate:.2f} - {rate_max:.2f}, Температура: {max_temp - 273.15:.2f}°C, Время: {t_time:.2f}c, Время между выстрелами {interval_between_shots:.2f}")
        if max_temp < T_interes_K:
            rate_min = rate
        else:
            rate_max = rate
            
        last_rate = rate
    
    # cur_temperatures = cur_temperatures_to_return
    if solution_found:
        dict1["str"] = f"Достижимо при \nскорострельности {round(last_rate)} выстр/мин"
        dict1["status"] = "success"
    elif ct == 1 and temperature_exceeded:
        dict1["str"] = "Достижимо при любой\nскорострельности"
        dict1["status"] = "success"
    else:
        dict1["str"] = f"Недостижимо при таких условиях.\nПоследняя T(r1) = {(max_temp - 273.15):.2f}°C\nПри скорострельности\n{round(last_rate)} выстр/мин\nВремя процесса {t_time:.2f}с"
        dict1["status"] = "failed"
    
    dict1["last_T"] = cur_temperatures
    dict1["all_T"] = all_T
    dict1["radii_nodes"] = radii_nodes
    dict1["ttime"] = ttime[:-1] if ttime else []
    dict1["final_rate"] = round(last_rate)
    dict1["final_temp"] = cur_temperatures[0] - 273.15
    dict1["final_time"] = t_time
    
    print("D")
    
    return dict1

def validate(v1, v2, v3, v4):
    if not((v1 > T_air_K) and (v2 > T_air_C) and (5000 <= v3 <= 25000) and (0 < v4 < 1200)):
        raise ValueError

if __name__ == '__main__':
    def on_viewport_resize():
        viewport_height = dpg.get_viewport_height()
        viewport_width = dpg.get_viewport_width()
        
        if viewport_width < MIN_WIDTH or viewport_height < MIN_HEIGHT:
            dpg.set_viewport_width(MIN_WIDTH)
            dpg.set_viewport_height(MIN_HEIGHT)
            return
        
        # Высота заголовка окна + отступы
        header_and_padding = 80
        # Доступная высота для содержимого главного окна
        available_height = viewport_height - header_and_padding
        plot_height = available_height // 2 - 5  # -5 для небольшого зазора
        if plot1_id:
            dpg.set_item_height(plot1_id, plot_height)
        if plot2_id:
            dpg.set_item_height(plot2_id, plot_height)

    # Глобальные ID графиков
    plot1_id = None
    plot2_id = None

    # Инициализация Dear PyGui
    dpg.create_context()
    
    # Шрифт
    with dpg.font_registry() :
        with dpg.font("C:/Windows/Fonts/arial.ttf", 16, tag="def_font"):
            dpg.add_font_range_hint(dpg.mvFontRangeHint_Cyrillic)
    
    # Создаём окно
    with dpg.window(label="Main Window", width=MIN_WIDTH, height=MIN_HEIGHT) as main_window:
        # Горизонтальный группировщик для левой и правой части
        with dpg.group(horizontal=True):
            # Лево
            with dpg.child_window(width=250, height=-1, no_scrollbar=True):
                dpg.add_text("T газов (К)")
                input1 = dpg.add_input_text(label="", default_value="1075", width=235)
                dpg.add_text("T(r1) желаемая (°C)")
                input2 = dpg.add_input_text(label="", default_value="220", width=235)
                dpg.add_text("Коэф. теплоотдачи газа (Вт/(м^2*К)")
                input3 = dpg.add_input_text(label="", default_value="19600", width=235)
                dpg.add_text("Время расчёта (с)")
                input4 = dpg.add_input_text(label="", default_value="300", width=235)

                dpg.add_spacer(height=10)
                output_text = dpg.add_text("Результат: (нажмите кнопку)")

                def on_button_click():
                    try:
                        dpg.disable_item("bt")
                        dpg.configure_item("bt", label="В процессе...")
                        v1 = float(dpg.get_value(input1))
                        v2 = float(dpg.get_value(input2))
                        v3 = float(dpg.get_value(input3))
                        v4 = float(dpg.get_value(input4))
                        validate(v1, v2, v3, v4)
                        dict1 = calc(v1, v2, v3, v4)
                        dpg.set_value(output_text, dict1["str"])
                        
                        new_x2 = [i*1000 for i in dict1["radii_nodes"]]
                        new_y2 = [i - 273.15 for i in dict1["last_T"]]
                        dpg.configure_item("series2", x=new_x2, y=new_y2)
                        dpg.fit_axis_data("x_axis2")
                        dpg.fit_axis_data("y_axis2")
                        
                        new_y1 = [i - 273.15 for i in dict1["all_T"]]
                        new_x1 = [i for i in dict1["ttime"]]
                        dpg.configure_item("series1", x=new_x1, y=new_y1)
                        dpg.fit_axis_data("x_axis1")
                        dpg.fit_axis_data("y_axis1")
                    
                    except ValueError:
                        dpg.set_value(output_text, "Ошибка: введите числа!\nИ проверьте их адекватность.")
                    
                    finally:
                        dpg.enable_item("bt")
                        dpg.configure_item("bt", label="Вычислить")

                dpg.add_button(label="Вычислить", callback=on_button_click, height=60, width=235, tag="bt")

            # Право
            with dpg.child_window(width=-1, height=-1, no_scrollbar=True):
                # Первый график
                plot1_id = dpg.generate_uuid()
                with dpg.plot(label="Температура на внутренней стенке ствола", tag=plot1_id, width=-1):
                    x_axis1 = dpg.add_plot_axis(dpg.mvXAxis, label="Время, c", tag="x_axis1")
                    y_axis1 = dpg.add_plot_axis(dpg.mvYAxis, label="Температура, °C", tag="y_axis1")
                    dpg.add_line_series(x=[], y=[], parent=y_axis1, tag="series1")

                dpg.add_spacer(height=10)

                # Второй график
                plot2_id = dpg.generate_uuid()
                with dpg.plot(label="Распределение T по радиусу ствола в последний момент времени", tag=plot2_id, width=-1):
                    x_axis2 = dpg.add_plot_axis(dpg.mvXAxis, label="Радиус, мм", tag="x_axis2")
                    y_axis2 = dpg.add_plot_axis(dpg.mvYAxis, label="Температура, °C", tag="y_axis2")
                    dpg.add_line_series(x=[], y=[], parent=y_axis2, tag="series2")
            
    dpg.bind_font("def_font")
    
    with dpg.theme() as global_theme:
        with dpg.theme_component(dpg.mvAll):
            dpg.add_theme_style(dpg.mvPlotStyleVar_LineWeight, 3.0, category=dpg.mvThemeCat_Plots)
            dpg.add_theme_style(dpg.mvStyleVar_FrameBorderSize, 1.0, category=dpg.mvThemeCat_Core)
    dpg.bind_theme(global_theme)
    # dpg.show_style_editor()

    dpg.create_viewport(title='Heat transfer modeling', width=1050, height=650, small_icon="pu.ico", large_icon="pu.ico")
    dpg.setup_dearpygui()
    dpg.set_viewport_resize_callback(on_viewport_resize)
    dpg.show_viewport()
    dpg.set_primary_window(main_window, True)
    dpg.start_dearpygui()
    dpg.destroy_context()
