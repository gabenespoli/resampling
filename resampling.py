import pandas as pd

import utils

filename = "./TutorialData.csv"
value_col = "Rating"
n_max = 100
iterations = 500
random_state = None
ci_method = "normal_range"
verbose = False

value_col_centered = f"{value_col}_centered"

df = pd.read_csv(filename)
df[value_col_centered] = df[value_col] - df[value_col].mean()

data = list()

for iter in range(1, iterations + 1):
    if verbose:
        print(f"Iteration {iter}")

    for n in range(1, n_max + 1):
        if verbose:
            print(".", end="")
        samples = df[value_col_centered].sample(
            n=n,
            replace=True,
            random_state=random_state,
        )
        mean_n = samples.mean()
        ci_n = utils.get_ci(samples, method=ci_method)
        data.append(
            {
                "iter": iter,
                "n": n,
                "mean": mean_n,
                "ci_lower": ci_n[0],
                "ci_upper": ci_n[1],
            }
        )

df_resampled = pd.DataFrame(data)

fig = utils.plot(df_resampled)
fig.show()
