import matplotlib.pyplot as plt
from tensorboard.backend.event_processing import event_accumulator


def plot_tensorboard_timeseries(event_file, tag, epochs=120):
    # Load TensorBoard event file
    ea = event_accumulator.EventAccumulator(event_file)
    ea.Reload()

    # Extract the time series data
    if tag not in ea.Tags()["scalars"]:
        raise ValueError(f"Tag '{tag}' not found in the event file.")

    events = ea.Scalars(tag)
    steps = [e.step / 54 for e in events][:epochs]
    values = [e.value for e in events][:epochs]
    total_time = events[epochs - 1].wall_time - events[0].wall_time
    return steps, values, total_time


def plot_scores(stepsLists, valueLists, wall_times, models):

    # Plot the time series
    plt.figure(figsize=(10, 5))
    for model_steps, value_values, wall_time, model in zip(
        stepsLists, valueLists, wall_times, models
    ):
        print(
            f"{model} wall time was {wall_time:.0f} seconds ({wall_time/3600:.2f} hours)"
        )
        # print(model_steps, value_values, model)
        plt.plot(model_steps, value_values, label=model)
    plt.xlabel("Epoch")
    plt.ylabel("Area under PR-curve")
    plt.legend()
    plt.grid(True)
    plt.savefig("au_pr_score.png")


# Usage example
event_file = "event_files/events.out.tfevents.1730291782.gcn30.local.snellius.surf.nl.2180618.0"  # Update with your TensorBoard event file path
event_file_dict = {
    "DinoV2": "event_files/events.out.tfevents.1730291782.gcn30.local.snellius.surf.nl.2180618.0",
    "EfficientNet": "event_files/events.out.tfevents.1730302269.gcn30.local.snellius.surf.nl.2222324.0",
}

# Extract data for each model
stepsLists, valuesLists, wall_times, models = zip(
    *[
        (
            *plot_tensorboard_timeseries(
                event_file_dict[model], "au_pr_curve/validation", epochs=120
            ),
            model,
        )
        for model in event_file_dict.keys()
    ]
)


plot_scores(stepsLists, valuesLists, wall_times, models)
