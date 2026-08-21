# Interactive Normative Charts

Vertex-resolution normative trajectories from the pretrained SNM-1000 model. The
visualization renders median cortical thickness across the cortical surface at a
selectable age, alongside lifespan charts for three exemplary vertices in insular,
prefrontal, and occipital cortex.

<video id="viz-demo" controls autoplay muted loop playsinline preload="metadata"
       style="display:block; margin:0 auto; max-width:500px; width:100%;
              border:1px solid var(--md-default-fg-color--lightest); border-radius:4px;">
  <source src="/normative_brain_charts/assets/visualization.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>

<script>
  document.querySelectorAll('#viz-demo').forEach(v => { v.playbackRate = 2.0; });
</script>

[**Open the interactive visualization**](/normative_brain_charts/code/notebooks/10_interactive_visualization/08_01_interactive_visualization.html){target=_blank}

!!! note
    The visualization is a self-contained HTML file of approximately 50 MB and may
    take a short while to load.