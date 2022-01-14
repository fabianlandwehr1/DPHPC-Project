const {plot, presets} = require("node-plotly.js");
const _ = require("lodash");

const colors = [
    '#1f77b4',  // muted blue
    '#ff7f0e',  // safety orange
    '#2ca02c',  // cooked asparagus green
    '#d62728',  // brick red
    '#9467bd',  // muted purple
    '#8c564b',  // chestnut brown
    '#e377c2',  // raspberry yogurt pink
    '#7f7f7f',  // middle gray
    '#bcbd22',  // curry yellow-green
    '#17becf'   // blue-teal
];

const symbols = [
    "circle",
    "cross",
    "star",
    "square",
    "diamond",
];


// Input data
const kernels = ["azimint", "conv2d", "durbin", "cavtflow", "gramschm"];
const I = {     // flops/byte
    "azimint": [499.875812094, 499.875812094, 374.9382848576],
    "conv2d": [1, 2, 231.36206],
    "durbin": [4001.00, 4001.00],
    "gramschm": [0.9965, 1.98312],
    "cavtflow": [1.71, 3.398],
}
const P = {     // flops/cycle
    "azimint": [1.3332029269, 3.9419231215, 1315.4345317308],
    "conv2d": [1.46, 1.8, 10.25],
    "durbin": [0.3214411657, 0.3331459966],
    "gramschm": [0.19, 14.665],
    "cavtflow": [6.269, 6.48],
}
const beta = 4 * 21.33 * 1024 * 1024 * 1024 / 200000000;    // bytes/cycle
const pi = {    // flops/cycle
    "azimint": 1315.4345317308 * 1.47,
    "conv2d": 348.977,
    "durbin": 0.3214411657 * 81,
    "gramschm": 1137,
    "cavtflow": 194.35,
}

const rangeLogX = [-2, 4];
const rangeLogY = [-1, 3.48];
const rangeX = rangeLogX.map(v => Math.pow(10, v));
const rangeY = rangeLogY.map(v => Math.pow(10, v));

// Memory bound line
const memBound = {
    x: rangeX,
    y: rangeX.map(x => beta * x),
    mode: "lines",
    line: {color: "#7f7f7f"},
};

// Memory bound label
const memBoundLabel = {
    x: Math.log10(0.0475),
    y: Math.log10(0.059 * beta),
    xanchor: "left",
    yanchor: "bottom",
    text: "memory bound",
    textangle: -38.5,
    font: {
        color: "#7f7f7f",
    },
    showarrow: false,
};

// Compute bound lines
const compBounds = kernels.map((k, i) => ({
    x: rangeX,
    y: [pi[k], pi[k]],
    mode: "lines",
    line: {
        color: colors[i],
    },
}));

// Compute bound labels
const yanchorBound = {
    "cavtflow": "top",
    "conv2d": "top",
};
const compBoundLabels = kernels.map((k, i) => ({
    x: rangeLogX[1],
    y: Math.log10(pi[k]),
    xanchor: "right",
    yanchor: yanchorBound[k] || "bottom",
    text: k,
    font: {
        color: colors[i],
    },
    showarrow: false,
}));

// Markers for a given kernel
function markers(name, index) {
    return {
        x: I[name],
        y: P[name],
        mode: "markers",
        marker: {
            color: colors[index],
            symbol: symbols[index],
            size: 8,
        },
    };
}

// Labels for markers
const xanchor = {
    // "azimint/opt2": "right",
    "cavtflow/base": "right",
    "conv2d/base": "right",
    "durbin/base": "right",
    "durbin/opt1": "right",
};
const yanchor = {
    "azimint/opt2": "bottom",
    "durbin/opt1": "bottom",
};
function labels(name, index) {
    return I[name].map((intensity, i) => {
        const text = `${ name }/${ !i ? "base" : "opt" + i }`;
        return {
            x: Math.log10(intensity),
            y: Math.log10(P[name][i]),
            xanchor: xanchor[text] || "left",
            yanchor: yanchor[text] || "top",
            text: text,
            font: {
                color: colors[index],
            },
            showarrow: false,
        };
    })
}

// Values of ticks
function tickvals(range) {
    return _.flatten(_.range(...range).map(start =>
        _.range(Math.pow(10, start), Math.pow(10, start + 1), Math.pow(10, start))));
}

// Tick labels
function ticktext(range) {
    return _.flatten(_.range(...range).map(start => [`10<sup>${start}</sup>`, ...Array(8).fill("")]));
}

plot([memBound, ...compBounds, ...kernels.map(markers)], _.defaultsDeep({
    width: 800,
    height: 560,
    title: {
        text: "Roofline Model - Xilinx Virtex UltraScale+ VU9P<br />Performance [flops/cycle]",
    },
    xaxis: {
        type: "log",
        range: rangeLogX,
        tickvals: [...tickvals(rangeLogX), 10000],
        ticktext: [...ticktext(rangeLogX), "10<sup>4</sup>"],
        tickangle: 0,
        gridwidth: 1,
        showgrid: true,
        zeroline: true,
        gridcolor: "#F0F0F0",
        title: {
            text: "Operational Intensity [flops/byte]",
        },
    },
    yaxis: {
        type: "log",
        range: rangeLogY,
        tickvals: tickvals(rangeLogY),
        ticktext: ticktext(rangeLogY),
        gridcolor: "#F0F0F0",
    },
    annotations: [memBoundLabel, ...compBoundLabels, ..._.flatten(kernels.map(labels))],
}, presets.asl));