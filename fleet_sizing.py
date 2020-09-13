import numpy as np
import pandas as pd
import json
import ast
import copy

class FleetSizing:

    # Parameters
    f = 1.5
    capacity = 4
    start_hour = 6
    frac = 1.0
    seed = 12345

    output_file = "D:/Dropbox (GaTech)/LeapHi/Pandemic routing/ODMTS designs/shuttlemile_1_0_transf_3_alpha_0_05_multipl_5_0_busfreqs_2_3_4_output.json"
    durations_file = "D:/Dropbox (Persoonlijk)/Postdoc/Code/AnalysisODMTS/data/duration_matrix_jan1_oct31.csv"
    trip_samples_file = "D:/Dropbox (Persoonlijk)/Postdoc/Code/AnalysisODMTS/data/sample_n_33000_odx.csv"
    locations_file = "D:/Dropbox (Persoonlijk)/Postdoc/Code/AnalysisODMTS/data/clustered_stops_1500feet_jan1_oct31.csv"

    # Output
    output = None
    shuttle_legs = None
    tasks = None
    shuttles = None
    hubs = None
    locations = None
    start_timestamp = None
    end_timestamp = None
    durations = None

    def __init__(self):
        pass

    def run(self):


        ################
        # PREPARE DATA #
        ################

        # Import ODMTS output file
        with open(self.output_file) as file:
            output = json.load(file)


        # Import duration matrix
        durations = pd.read_csv(self.durations_file)
        durations['stop_id'] = durations['stop_id'].astype(int).astype(str)
        durations.set_index('stop_id', inplace=True)
        durations.columns = durations.columns.astype(float).astype(int).astype(str)
        assert all(durations.index == durations.columns)
        durations = durations.applymap(lambda x: int(np.ceil(60 * x))) # round up to assure triangle inequality

        # Get hubs and locations
        hubs = pd.DataFrame(output['design']['hubs'].keys(), columns=['stop_id'])

        locations = pd.read_csv(self.locations_file,
                                dtype={'stop_id': 'int',
                                       'stop_lat': 'float',
                                       'stop_lon': 'float'})

        locations['stop_id'] = locations['stop_id'].astype(str)
        locations.set_index('stop_id', inplace=True)


        # Use trip splitting to replace trips by sampled trips
        samples = pd.read_csv(self.trip_samples_file,
                             dtype={'start_cluster': str, 'end_cluster': str, 'count': int, 'start_times': str},
                             index_col=0).reset_index(drop=True)
        samples['start_times'] = samples['start_times'].apply(ast.literal_eval)
        samples = samples.sample(frac=self.frac, random_state=self.seed)

        trips = list()

        for index, sample in samples.iterrows():

            trip = dict()
            trip['passengers'] = int(sample['count'])
            trip['origin_stop_id'] = sample['start_cluster']
            trip['destination_stop_id'] = sample['end_cluster']
            trip['departure_times'] = sample['start_times']

            trip['legs'] = list()
            for leg in output['trip_splitting'][trip['origin_stop_id']][trip['destination_stop_id']]['legs']:
                if 'i' in leg:
                    trip['legs'].append(output['design']['legs'][leg['i']].copy())
                else:
                    trip['legs'].append(
                                {'board_stop_id': leg['b'],
                                 'alight_stop_id': leg['a'],
                                 'mode': 0,
                                 'travel_time': durations.loc[leg['b'], leg['a']],
                                 'total_time': durations.loc[leg['b'], leg['a']]}
                            )

            trips.append(trip)

        output['trips'] = trips


        # Determine start_timestamp and initialize end_timestamp
        start_timestamp = pd.to_datetime(output['trips'][0]['departure_times'][0])
        start_timestamp = start_timestamp.replace(hour=self.start_hour, minute=0, second=0)
        start_timestamp = int(start_timestamp.timestamp()) - 3600  # start fleet sizing earlier to allow early pickups
        end_timestamp = start_timestamp


        # add start_times and end_times to trip legs
        for trip in output['trips']:
            times = [int(i.timestamp()) - start_timestamp for i in pd.to_datetime(trip['departure_times'])]

            for leg in trip['legs']:
                leg['start_times'] = times.copy()
                times = [i + leg['total_time'] for i in times]
                leg['end_times'] = times.copy()

                for end_time in leg['end_times']:
                    if start_timestamp + end_time + durations.loc[leg['board_stop_id'], leg['alight_stop_id']] > end_timestamp:
                        end_timestamp = start_timestamp + end_time + durations.loc[leg['board_stop_id'], leg['alight_stop_id']]


        # Prepare shuttle_Legs
        shuttle_legs = pd.DataFrame([leg for trip in output['trips'] for leg in trip['legs'] if leg['mode'] == 0])
        shuttle_legs.loc[:, 'board_at_hub'] = shuttle_legs.loc[:, 'board_stop_id'].apply(
            lambda x: bool((hubs == x).sum()['stop_id']))
        shuttle_legs.loc[:, 'alight_at_hub'] = shuttle_legs.loc[:, 'alight_stop_id'].apply(
            lambda x: bool((hubs == x).sum()['stop_id']))

        shuttle_legs = shuttle_legs.explode(column='start_times').reset_index(drop=True)
        shuttle_legs.rename(columns={'start_times': 'request_open'}, inplace=True)

        def assign_to_hub(row):
            if row['board_at_hub']:
                return row['board_stop_id']
            elif row['alight_at_hub']:
                return row['alight_stop_id']
            else:
                return hubs.loc[hubs['stop_id'].apply(lambda x: durations[x][row['board_stop_id']]).idxmin(), 'stop_id']

        shuttle_legs.loc[:, 'hub_assignment'] = shuttle_legs.apply(assign_to_hub, axis=1)
        shuttle_legs.loc[:, 'non_hub'] = shuttle_legs.apply(lambda row:
                                                            row['board_stop_id'] if row['hub_assignment'] == row[
                                                                'alight_stop_id'] else row['alight_stop_id'],
                                                            axis=1)
        shuttle_legs.loc[:, 'type'] = shuttle_legs.apply(lambda row:
                                                         'pickup' if row['hub_assignment'] == row['alight_stop_id'] else
                                                         'dropoff' if row['hub_assignment'] == row['board_stop_id'] else
                                                         'direct', axis=1)
        shuttle_legs = shuttle_legs[['hub_assignment', 'non_hub', 'type', 'request_open', 'travel_time']]


        # Print simple statistics and remove direct trips
        direct_legs = shuttle_legs[shuttle_legs['type'] == 'direct']
        hub_legs = shuttle_legs[(shuttle_legs['type'] == 'pickup') | (shuttle_legs['type'] == 'dropoff')]

        print("Direct legs: " + str(len(direct_legs)) + ", total travel time: " + str(direct_legs['travel_time'].sum()) + " s")
        print("Hub legs: " + str(len(hub_legs)) + ", total travel time: " + str(hub_legs['travel_time'].sum()) + " s")
        print("Scale factor hub only -> all (#legs based): " + str(len(shuttle_legs)/len(hub_legs)))
        print("Scale factor hub only -> all (time based): " + str(shuttle_legs['travel_time'].sum()/hub_legs['travel_time'].sum()))

        shuttle_legs = shuttle_legs[(shuttle_legs['type'] == 'pickup') | (shuttle_legs['type'] == 'dropoff')]

        counts = shuttle_legs.groupby('hub_assignment').count()
        min_hub = counts.idxmin()[0]
        max_hub = counts.idxmax()[0]
        print("Smallest group: " + min_hub + " -> " + str(counts.loc[min_hub,:][0]))
        print("Average group: " + str(counts.mean()[0]))
        print("Largest group: " + max_hub + " -> " + str(counts.loc[max_hub,:][0]))



        #################################
        # GREEDY CONSTRUCTION HEURISTIC #
        #################################

        global_lb = 0
        global_ub = end_timestamp - start_timestamp


        # Prepare tasks based on the shuttle legs
        tasks = shuttle_legs.copy()

        def create_task(row):
            travel_time_from_hub = durations.loc[row['hub_assignment'], row['non_hub']]
            travel_time_to_hub = durations.loc[row['non_hub'], row['hub_assignment']]

            ftt_from = int(np.round(self.f * travel_time_from_hub))
            ftt_to = int(np.round(self.f * travel_time_to_hub))

            if row['type'] == 'pickup':
                start_hub_tw = (global_lb,
                                row['request_open'] + ftt_to - travel_time_to_hub - travel_time_from_hub)
                non_hub_tw = (row['request_open'],
                              row['request_open'] + ftt_to - travel_time_to_hub)
                end_hub_tw = (row['request_open'] + travel_time_to_hub,
                              row['request_open'] + ftt_to)

            elif row['type'] == 'dropoff':
                start_hub_tw = (row['request_open'],
                                row['request_open'] + ftt_from - travel_time_from_hub)
                non_hub_tw = (row['request_open'] + travel_time_from_hub,
                              row['request_open'] + ftt_from)
                end_hub_tw = (row['request_open'] + travel_time_from_hub + travel_time_to_hub,
                              global_ub)

            return start_hub_tw[0], start_hub_tw[1], non_hub_tw[0], non_hub_tw[1], end_hub_tw[0], end_hub_tw[1]

        tasks[['start_hub_lb', 'start_hub_ub', 'non_hub_lb', 'non_hub_ub', 'end_hub_lb', 'end_hub_ub']] = \
            tasks.apply(create_task, axis=1, result_type='expand')

        tasks = tasks[
            ['hub_assignment', 'non_hub', 'type', 'start_hub_lb', 'start_hub_ub', 'non_hub_lb', 'non_hub_ub', 'end_hub_lb',
             'end_hub_ub']]


        # Sort tasks by the deadline for the vehicle leaving the depot
        tasks.sort_values('start_hub_ub', inplace=True)


        # Define helper functions
        def reduce_time_windows(path, task_types, time_windows, capacity, check_capacity=True):

            # check capacity constraint (only correct if path is a single route from depot to depot, otherwise load overestimated)
            if check_capacity:
                load = sum(1 for task_type in task_types if task_type == 'dropoff')
                for i in range(0, len(task_types)):
                    if task_types[i] == 'dropoff':
                        load -= 1
                    elif task_types[i] == 'pickup':
                        load += 1
                    else:
                        assert (pd.isnull(task_types[i]))
                    if load > capacity:
                        return False
                    assert (load >= 0)

            # reduce time windows
            for i in range(1, len(time_windows)):
                time_windows[i][0] = max(time_windows[i - 1][0] + durations.loc[path[i - 1], path[i]], time_windows[i][0])
            for i in reversed(range(0, len(time_windows) - 1)):
                time_windows[i][1] = min(time_windows[i + 1][1] - durations.loc[path[i], path[i + 1]], time_windows[i][1])

            # check time window constraints
            for i in range(0, len(time_windows)):
                if time_windows[i][0] > time_windows[i][1]:
                    return False

            return True


        def get_best_insert_position(task_insert, path, task_types, time_windows, capacity):

            best_position = None
            best_end_hub_lb = None

            for position in range(1, len(path)):

                new_path = path[:position] + [task_insert.non_hub] + path[position:]
                new_task_types = task_types[:position] + [task_insert.type] + task_types[position:]
                new_time_windows = [[max(time_windows[0][0], task_insert.start_hub_lb),
                                     min(time_windows[0][1], task_insert.start_hub_ub)]] + \
                                   copy.deepcopy(time_windows[1:position]) + \
                                   [[task_insert.non_hub_lb, task_insert.non_hub_ub]] + \
                                   copy.deepcopy(time_windows[position:-1]) + \
                                   [[max(time_windows[-1][0], task_insert.end_hub_lb),
                                     min(time_windows[-1][1], task_insert.end_hub_ub)]]

                feasible = reduce_time_windows(new_path, new_task_types, new_time_windows, capacity)

                if feasible:
                    end_hub_lb = new_time_windows[-1][0]
                    if best_end_hub_lb == None or end_hub_lb < best_end_hub_lb:
                        best_position = position
                        best_end_hub_lb = end_hub_lb

            if best_position == None:
                return None
            else:
                return (best_position, best_end_hub_lb)


        def insert_task_at_position(task_insert, position, path, task_types, task_indices, time_windows, nb_pickups, nb_dropoffs, capacity):

            new_path = path[:position] + [task_insert.non_hub] + path[position:]
            new_task_types = task_types[:position] + [task_insert.type] + task_types[position:]
            new_task_indices = task_indices[:position] + [task_insert.Index] + task_indices[position:]
            new_time_windows = [[max(time_windows[0][0], task_insert.start_hub_lb),
                                 min(time_windows[0][1], task_insert.start_hub_ub)]] + \
                               copy.deepcopy(time_windows[1:position]) + \
                               [[task_insert.non_hub_lb, task_insert.non_hub_ub]] + \
                               copy.deepcopy(time_windows[position:-1]) + \
                               [[max(time_windows[-1][0], task_insert.end_hub_lb),
                                 min(time_windows[-1][1], task_insert.end_hub_ub)]]

            new_nb_pickups = nb_pickups
            new_nb_dropoffs = nb_dropoffs
            if task_insert.type == 'pickup':
                new_nb_pickups += 1
            else:
                new_nb_dropoffs += 1

            feasible = reduce_time_windows(new_path, new_task_types, new_time_windows, capacity)
            assert(feasible)

            return new_path, new_task_types, new_task_indices, new_time_windows, new_nb_pickups, new_nb_dropoffs


        # Perform algorithm

        #keep track of all shuttles
        shuttles = pd.DataFrame(columns=['hub_assignment', 'time', 'path', 'task_types', 'task_indices'])

        # Progress bar for use in Jupyter notebooks
        # progress = 0
        # progress_bar = tqdm.notebook.tqdm(total=len(tasks))

        for hub_assignment in tasks['hub_assignment'].unique():

            #keep local copies of tasks and shuttles for this specific hub
            tasks_hub = tasks[(tasks['hub_assignment'] == hub_assignment)].copy()
            shuttles_hub = pd.DataFrame(pd.Series(
                {'hub_assignment': hub_assignment, 'time': global_lb, 'path': list(), 'task_types': list(), 'task_indices': list()}
                )).transpose()
            shuttles_hub['time'] = shuttles_hub['time'].astype(int)

            #keep set which tasks are already scheduled
            scheduled = set()

            #keep scheduling the first next task, and add as many other tasks as possible
            for task_index, task in tasks_hub.iterrows():

                #stop when everything is scheduled
                if len(scheduled) == len(tasks_hub):
                    break

                #skip tasks that have been scheduled before
                if task_index in scheduled:
                    continue

                # progress_bar.n = progress + len(scheduled)
                # progress_bar.refresh()

                #select next available shuttle
                if len(shuttles_hub) == 1:
                    shuttle = shuttles_hub.loc[0, :].copy()
                else:
                    shuttle = shuttles_hub.loc[shuttles_hub['time'].idxmin(), :].copy()

                #infeasible, new shuttle needed
                if shuttle['time'] > task['start_hub_ub']:
                    shuttles_hub = shuttles_hub.append(
                        {'hub_assignment': hub_assignment, 'time': global_lb, 'path': list(), 'task_types': list(),
                         'task_indices': list()}, ignore_index=True)
                    shuttle = shuttles_hub.iloc[-1, :].copy()

                # advance to earliest departure time
                shuttle['time'] = max(task['start_hub_lb'], shuttle['time'])

                # initialize route
                path = [hub_assignment, task['non_hub'], hub_assignment]
                task_types = [np.nan, task['type'], np.nan]
                task_indices = [np.nan, task_index, np.nan]
                time_windows = [[shuttle['time'], task['start_hub_ub']],
                                [task['non_hub_lb'], task['non_hub_ub']],
                                [task['end_hub_lb'], task['end_hub_ub']]]
                reduce_time_windows(path, task_types, time_windows, self.capacity)
                nb_pickups = 0
                nb_dropoffs = 0
                if task['type'] == 'pickup':
                    nb_pickups += 1
                else:
                    nb_dropoffs += 1

                scheduled.add(task_index)

                # Keep adding tasks until the route is full
                while True:

                    best_task_insert = None
                    best_position = None
                    best_end_hub_lb = None

                    # Insert a single task
                    for task_insert in tasks_hub.itertuples():

                        if task_insert.Index in scheduled:
                            continue

                        if best_end_hub_lb != None and task_insert.end_hub_lb > best_end_hub_lb:
                            continue

                        if nb_pickups == self.capacity and task_insert.type == 'pickup':
                            continue

                        if nb_dropoffs == self.capacity and task_insert.type == 'dropoff':
                            continue

                        if task_insert.end_hub_ub < time_windows[-1][0] or task_insert.end_hub_lb > time_windows[-1][1]:
                            continue

                        if task_insert.start_hub_ub < time_windows[0][0] or task_insert.start_hub_lb > time_windows[0][1]:
                            continue

                        result = get_best_insert_position(task_insert, path, task_types, time_windows, self.capacity)

                        if result == None:
                            continue
                        else:
                            position = result[0]
                            end_hub_lb = result[1]
                            if best_end_hub_lb == None or end_hub_lb < best_end_hub_lb:
                                best_task_insert = task_insert
                                best_position = position
                                best_end_hub_lb = end_hub_lb

                    # No additions found, break
                    if best_position == None:
                        break

                    # Add new task and update data
                    path, task_types, task_indices, time_windows, nb_pickups, nb_dropoffs = \
                        insert_task_at_position(best_task_insert, best_position, path, task_types, task_indices, time_windows,
                                                nb_pickups, nb_dropoffs, self.capacity)

                    # Register that the task is scheduled
                    scheduled.add(best_task_insert.Index)

                    # Also stop if the route is now full
                    if nb_pickups == self.capacity and nb_dropoffs == self.capacity:
                        break

                shuttles_hub.at[shuttle.name, 'time'] = int(time_windows[-1][0])
                shuttles_hub.at[shuttle.name, 'path'] = shuttles_hub.loc[shuttle.name, 'path'] + path
                shuttles_hub.at[shuttle.name, 'task_types'] = shuttles_hub.loc[shuttle.name, 'task_types'] + task_types
                shuttles_hub.at[shuttle.name, 'task_indices'] = shuttles_hub.loc[shuttle.name, 'task_indices'] + task_indices

            # progress += len(tasks_hub)

            #add local shuttles_hub to all shuttles
            shuttles = pd.concat([shuttles, shuttles_hub]).reset_index(drop=True)

        # progress_bar.close()


        # Validation of the shuttle routes
        def get_departure_times(row):

            time_windows = [[global_lb, global_ub] for i in range(len(row['path']))]
            hub_visits = [i for i in range(len(row['task_indices'])) if pd.isnull(row['task_indices'][i])]

            for i in range(1, len(time_windows) - 1):
                if not pd.isnull(row['task_indices'][i]):
                    pr = max(j for j in hub_visits if j < i)
                    ne = min(j for j in hub_visits if j > i)

                    time_windows[pr][0] = max(time_windows[pr][0], tasks.loc[row['task_indices'][i], 'start_hub_lb'])
                    time_windows[pr][1] = min(time_windows[pr][1], tasks.loc[row['task_indices'][i], 'start_hub_ub'])
                    time_windows[i][0] = max(time_windows[i][0], tasks.loc[row['task_indices'][i], 'non_hub_lb'])
                    time_windows[i][1] = min(time_windows[i][1], tasks.loc[row['task_indices'][i], 'non_hub_ub'])
                    time_windows[ne][0] = max(time_windows[ne][0], tasks.loc[row['task_indices'][i], 'end_hub_lb'])
                    time_windows[ne][1] = min(time_windows[ne][1], tasks.loc[row['task_indices'][i], 'end_hub_ub'])

            feasible = reduce_time_windows(row['path'], row['task_types'], time_windows, self.capacity, False)
            assert (feasible)

            time_windows[0][0] = time_windows[1][0] - durations.loc[row['path'][0], row['path'][1]]

            return [time_window[0] for time_window in time_windows]

        shuttles['departure_times'] = shuttles.apply(get_departure_times, axis=1)
        shuttles['arrival_times'] = shuttles.apply(lambda row:
                                                   [row['departure_times'][i - 1] + durations.loc[row['path'][i - 1], row['path'][i]]
                                                    if i > 0 else row['departure_times'][0] for i in range(len(row['path']))],
                                                   axis=1)

        def get_loads(row):

            task_indices = row['task_indices']
            n = len(task_indices)

            loads_on_arrival = [0 for i in range(n)]
            loads_on_departure = [0 for i in range(n)]

            hub_types = row['task_types'].copy()

            for i in range(n):
                if pd.isnull(task_indices[i]):
                    if i == 0 or pd.isnull(task_indices[i - 1]):
                        hub_types[i] = 1
                    else:
                        hub_types[i] = -1

            for i in range(n):
                if hub_types[i] == 1 or hub_types[i] == -1:
                    if i == 0 or hub_types[i] == 1:
                        j = [k for k in range(i + 1, n) if hub_types[k] == -1][0]
                        loads_on_departure[i] = sum([k == 'dropoff' for k in hub_types[(i + 1):j]])
                        loads_on_arrival[i] = 0
                    else:
                        j = [k for k in range(i) if hub_types[k] == 1][-1]
                        loads_on_arrival[i] = sum([k == 'pickup' for k in hub_types[j:i]])
                        loads_on_departure[i] = 0
                else:
                    if hub_types[i] == 'pickup':
                        loads_on_arrival[i] = loads_on_departure[i - 1]
                        loads_on_departure[i] = loads_on_arrival[i] + 1
                    else:
                        loads_on_arrival[i] = loads_on_departure[i - 1]
                        loads_on_departure[i] = loads_on_arrival[i] - 1

            return loads_on_arrival, loads_on_departure

        shuttles[['loads_on_arrival', 'loads_on_departure']] = shuttles.apply(get_loads, axis=1, result_type='expand')
        shuttles['nb_tasks'] = shuttles['task_indices'].apply(lambda x: (~pd.isnull(x)).sum())

        self.output = output
        self.shuttle_legs = shuttle_legs
        self.tasks = tasks
        self.shuttles = shuttles
        self.hubs = hubs
        self.locations = locations
        self.start_timestamp = start_timestamp
        self.end_timestamp = end_timestamp
        self.durations = durations