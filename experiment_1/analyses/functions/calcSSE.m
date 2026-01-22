function sse = calcSSE(measured_data, estimated_data)

    sse = sum((measured_data - estimated_data).^2);

end