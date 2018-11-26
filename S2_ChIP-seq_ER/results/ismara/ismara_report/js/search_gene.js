/*global $:false */

function myFormatter(label, series) {
    if (series.lines.show == 0) {
        return(null)
    } 
    return(label);
}

function get_expr_mean(data) {
    var sum = 0;
    for (var i in data["data"]) {
        sum = sum + parseFloat(data["data"][i][1])
    }
    mean = sum / data["data"].length
    return(mean);   
}


function get_min_max(data) {
    var min = data["data"][0][1];
    var max = data["data"][0][1];
    for (var i in data["data"]) {
        var value = data["data"][i][1];
        if (value > max) {
            max = value;
        } 
        else if (value < min) {
            min = value;
        }
    }
    return([min, max]);    
}

function make_expr_plot(data, samples, id, plots) {
    var factor = 1
    if (samples.length > 14) {
        factor = Math.floor(samples.length / 10);
    }
    var samples2show = [];
    for (i = 0; i < samples.length; i++) {
        if (i % factor == 0) {
            samples2show[i] = samples[i]
        } else {
            samples2show[i] = [i, ""];
        }
    }
    options = {legend: {container: $("#" + id.replace("_div", "_legend"))}, xaxis: {ticks: samples2show}};
    plots.push($.plot("#" + id, data, options));
    var i = 0;
    // fix colors in prom expr plot
    $.each(data, function(key, val) {
        val.color = i;
        ++i;
    });
}


function sum_datasets(data1, data2) {
    for (var i in data1["data"]) {
        data1["data"][i][1] = data1["data"][i][1] + data2["data"][i][1];
    }
    return(data1);
}

function subtract_datasets(data1, data2) {
    for (var i in data1["data"]) {
        data1["data"][i][1] = data1["data"][i][1] - data2["data"][i][1];
    }
    return(data1);
}

function changeStatus(event) {
    plotType = event.data.param1;
/*    var show = $(this).attr('show');
    if (show === "1" || show === undefined) {
        $(this).attr("show", "0");
        $(this).find("button").text("Off");
    }
    else {
        $(this).attr("show", "1");
        $(this).find("button").text("On");
    }*/
    if (plotType == "expr") {
        replot($(this));
    } else if (plotType == "act") {
        replot_act($(this));
    }
}

function replot(elem) {
    console.log(elem);
    console.log($(elem).text());
    var label = $(elem).text();
    var plot_idx = parseInt($(elem).parents("div").data("prom-num"));
    var plot_data = plots[plot_idx].getData();
    for (var i in plot_data) {
        if (plot_data[i].label == label) {
            plot_data[i].lines.show = !plot_data[i].lines.show;
            break;
        }
    }
    plots[plot_idx].setData(plot_data);
    plots[plot_idx].draw();
    return("ok");
}

function replot_act(elem) {
    var label = $(elem).data("mat-name");
    var state = $(elem).prop("checked");
    var plot_idx = parseInt($(elem).parents("table").data("prom-num"));
    var plot_data = plots[plot_idx].getData();
    var mat_expr_data;
    for (var i in plot_data) {
        if (plot_data[i].label == label) {
            mat_expr_data = plot_data[i];
            break;
        }
    }
    if (state) {
        plot_data[plot_data.length - 1] = sum_datasets(plot_data[plot_data.length - 1], mat_expr_data);
    } else {
        plot_data[plot_data.length - 1] = subtract_datasets(plot_data[plot_data.length - 1], mat_expr_data);
    }
    plot_data[plot_data.length - 1].lines.show = 1;
    plots[plot_idx].setData(plot_data);
    plots[plot_idx].draw();
    return("ok");
}
