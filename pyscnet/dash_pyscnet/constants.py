# Colors for contour
ternary_color = [
    "#FFDF7F",
    "#FFDF7F",
    "#FFDF7F",
    "#FFF9DC",
    "#FFF9DC",
    "#FFF9DC",
    "#FFF9DC",
    "#FFF9DC",
    "#FFF9DC",
    "#FFF2AD",
    "#FFF2AD",
    "#FFF2AD",
    "#FDFFE9",
    "#FDFFE9",
    "#FDFFE9",
    "#FDFFE9",
]

# Colors for legend
colors = [
    "#001f3f",
    "#0074d9",
    "#3d9970",
    "#111111",
    "#01ff70",
    "#ffdc00",
    "#ff851B",
    "#ff4136",
    "#85144b",
    "#f012be",
    "#b10dc9",
    "#AAAAAA",
    "#111111",
]

# Ternary axis contour
ternary_contour = {
    "Carbonate Dominated Lith": [
        (0, 100, 0),  # (Quartz, Carbonate, Clay)
        (20, 80, 0),
        (0, 80, 20),
    ],
    "Silica Dominated Lith": [
        (100, 0, 0),
        (80, 20, 0),
        (80, 0, 20),
    ],
    "Clay Dominated Lith": [
        (0, 20, 80),
        (0, 0, 100),
        (20, 0, 80),
    ],
    "Silica rich Carbonate Mudstone": [
        (20, 80, 0),
        (10, 80, 10),
        (40, 50, 10),
        (50, 50, 0),
    ],
    "Silica rich Argillaceous Mudstone": [
        (40, 10, 50),
        (10, 10, 80),
        (20, 0, 80),
        (50, 0, 50),
    ],
    "Carbonate-rich Argillaceous Mudstone": [
        (0, 50, 50),
        (0, 20, 80),
        (10, 10, 80),
        (10, 40, 50),
    ],
    "Carbonate rich Silliceous Mudstone": [
        (80, 10, 10),
        (50, 10, 40),
        (50, 0, 50),
        (80, 0, 20),
    ],
    "Clay rich Silliceous Mudstone": [
        (80, 20, 0),
        (80, 10, 10),
        (50, 40, 10),
        (50, 50, 0),
    ],
    "Clay-rich Carbonate Mudstone": [
        (10, 80, 10),
        (0, 80, 20),
        (0, 50, 50),
        (10, 50, 40),
    ],
    "Carbonate Silliceous Mudstone": [
        (50, 50, 0),
        (30, 50, 20),
        (50, 30, 20),
    ],
    "Argillaceous Carbonate Mudstone": [
        (20, 50, 30),
        (0, 50, 50),
        (20, 30, 50),
    ],
    "Argillaceous Sillicaous Mudstone": [
        (50, 20, 30),
        (30, 20, 50),
        (50, 0, 50),
    ],
    "Mixed Silliceous Mudstone": [
        (80, 10, 10),
        (50, 40, 10),
        (50, 10, 40),
    ],
    "Mixed Mudstone": [
        (30, 50, 20),
        (20, 50, 30),
        (20, 30, 50),
        (30, 20, 50),
        (50, 20, 30),
        (50, 30, 20),
    ],
}